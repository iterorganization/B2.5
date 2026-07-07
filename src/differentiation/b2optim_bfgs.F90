!=======================================================================
! b2optim_bfgs.f90
!
! Written by Claude Code
!
! Standalone bound-constrained BFGS optimiser for B2.5,
! structurally equivalent to b2optim_tao but with no PETSc/TAO
! dependency.  The embedded bfgs_bound_mod module provides the solver;
! everything else mirrors b2optim_tao as closely as possible:
!
!   - b2mn_step      is the cheap forward-only evaluator (line search)
!   - b2mn_step_diff   is the expensive forward+gradient evaluator
!                    (called once per accepted BFGS iterate)
!   - parameter rescaling, sigma/mean/shift/corr handling,
!     active-drift reset, intermediate state writing, and
!     TGT / ADJ gradient paths are all preserved
!
! Build (example):
!   gfortran -O2 [your b2mod objects/flags] -o b2optim_bfgs b2optim_bfgs.f90
!
!=======================================================================

!=======================================================================
! Dependency-free bound-constrained L-BFGS module
!=======================================================================
module bfgs_bound_mod
  implicit none
  private
  public :: func_if, func_grad_if, bfgs_bound_solve

  ! Cheap function-only evaluator (used inside the line search)
  abstract interface
     subroutine func_if(n, x, f)
       integer,          intent(in)  :: n
       double precision, intent(in)  :: x(n)
       double precision, intent(out) :: f
     end subroutine func_if
  end interface

  ! Expensive function+gradient evaluator (once per accepted iterate)
  abstract interface
     subroutine func_grad_if(n, x, f, g)
       integer,          intent(in)  :: n
       double precision, intent(in)  :: x(n)
       double precision, intent(out) :: f
       double precision, intent(out) :: g(n)
     end subroutine func_grad_if
  end interface

contains

  !---------------------------------------------------------------------
  ! bfgs_bound_solve
  !
  !  Minimizes f(x) s.t. l <= x <= u via bound-constrained L-BFGS with
  !  diagonal Broyden scaling of the initial Hessian (mirrors PETSc's
  !  MATLMVMBFGS default -mat_lmvm_scale_type diagonal, per Gilbert &
  !  Lemarechal 1989), rather than a single scalar gamma. This matters
  !  a lot when different parameters have very different curvature.
  !
  !  f_only   – cheap forward-only call  (line search trials)
  !  fg       – expensive forward+grad   (once per accepted step)
  !  memsize  – L-BFGS history length (TAO default = 1)
  !
  !  h0_type – which initial-Hessian (J0) scaling scheme to use:
  !          0 = SCALAR  (Nocedal gamma = (s^T y)/(y^T y) from the most
  !              recent pair; classic L-BFGS scaling, simple & robust)
  !          1 = DIAGONAL (per-component Broyden-class diagonal scaling,
  !              closer in spirit to PETSc's MATLMVMBFGS default, but a
  !              more elaborate recursive scheme -- see 'theta' below)
  !              Default is 0 (SCALAR).
  !
  !  Three convergence tolerances matching TaoSetTolerances (TAO fires
  !  when ANY is met; pass 1e30 to disable one). These match TAO's
  !  ACTUAL internal semantics (note grtol and gttol use different
  !  normalizations than their names might suggest):
  !   gatol  – absolute gradient norm:        ||g||        <= gatol
  !   grtol  – gradient relative to cost:     ||g||/|f|    <= grtol
  !   gttol  – gradient reduction from start: ||g||/||g0|| <= gttol
  !
  !  theta – (only used when h0_type=1) convex mixing factor for the
  !          diagonal J0 scaling, matching PETSc's
  !          -tao_bqnls_mat_lmvm_theta: Bdiag = (1-theta)*BFGS +
  !          theta*DFP diagonal updates. theta=0 is pure diagonal BFGS;
  !          PETSc's own default is 0.125 (found to be numerically
  !          fragile, use with caution).
  !
  !  alpha_init – initial trial step length for the backtracking line
  !          search at every iteration, matching TAO's -tao_ls_stepinit
  !          (TAO's own default is 1.0).
  !
  !  n_updates, n_rejects, n_resets – diagnostic counters mirroring
  !          TAO's "Mat Object" view (Number of updates/rejects/resets)
  !
  !  conv_reason – which criterion triggered convergence (only
  !          meaningful when info==0):
  !            0 = gatol   (||g|| <= gatol)
  !            1 = grtol   (||g||/|f| <= grtol)
  !            2 = gttol   (||g||/||g0|| <= gttol)
  !           -1 = did not converge (info /= 0)
  !  gnorm_final – the projected gradient norm ||g|| at exit, useful
  !          for reporting regardless of which criterion (or none)
  !          triggered the exit.
  !
  !  info:  0 = converged, 1 = max iterations, 2 = line search failure
  !---------------------------------------------------------------------
  subroutine bfgs_bound_solve(n, x, l, u, f_only, fg,          &
                               maxit, gatol, grtol, gttol,     &
                               memsize, h0_type, theta,        &
                               rescale_rho,                    &
                               alpha_init,                     &
                               iprint,                         &
                               iter, info, conv_reason,        &
                               gnorm_final,                    &
                               n_updates, n_rejects, n_resets)
    integer,          intent(in)    :: n, maxit, memsize, h0_type, iprint
    double precision, intent(inout) :: x(n)
    double precision, intent(in)    :: l(n), u(n)
    procedure(func_if)              :: f_only
    procedure(func_grad_if)         :: fg
    double precision, intent(in)    :: gatol, grtol, gttol, theta
    double precision, intent(in)    :: rescale_rho
    double precision, intent(in)    :: alpha_init
    integer,          intent(out)   :: iter, info, conv_reason
    double precision, intent(out)   :: gnorm_final
    integer,          intent(out)   :: n_updates, n_rejects, n_resets

    ! L-BFGS circular buffer
    double precision, allocatable :: S(:,:), Y(:,:), rho(:)
    double precision, allocatable :: g(:), gold(:), xold(:), d(:), q(:)
    double precision, allocatable :: Bdiag(:), H0diag(:)
    double precision, allocatable :: H0new(:), H0diag_stored(:)
    logical,          allocatable :: active(:), active_old(:)
    double precision, allocatable :: alpha_ls(:)

    double precision :: f, fold, alpha, beta_ls, gnorm, gnorm0
    double precision :: sty_free, yty_free, H0scale, sBs, pred, theta_c
    double precision :: bfgs_i, dfp_i, yHy, sigma_scale, rescale_rho_c
    double precision, parameter :: BTOL      = 1.0d-10
    double precision, parameter :: C1        = 1.0d-4
    double precision, parameter :: RHOLS     = 0.5d0
    double precision, parameter :: ALPHA_MIN = 1.0d-20
    double precision, parameter :: CURV_TOL  = 1.0d-10
    double precision, parameter :: BDIAG_MIN = 1.0d-8
    double precision, parameter :: BDIAG_MAX = 1.0d8
    double precision, parameter :: RESCALE_TOL = 1.0d-8  ! PETSc's ldb->tol
    integer :: i, j, k, m, ncur, newest, col
    logical :: active_set_changed

    m = max(1, memsize)
    theta_c = min(max(theta, 0.0d0), 1.0d0)   ! clamp to [0,1] as PETSc does
    rescale_rho_c = min(max(rescale_rho, 0.0d0), 1.0d0)
    allocate(S(n,m), Y(n,m), rho(m))
    allocate(g(n), gold(n), xold(n), d(n), q(n), alpha_ls(m))
    allocate(active(n), active_old(n), Bdiag(n), H0diag(n))
    allocate(H0new(n), H0diag_stored(n))
    S = 0.0d0;  Y = 0.0d0;  rho = 0.0d0
    Bdiag = 1.0d0           ! diagonal Hessian approx, starts at identity
    H0diag_stored = 1.0d0   ! corresponding inverse scaling, starts at identity
    ncur = 0   ! number of (s,y) pairs stored so far
    newest = 0 ! circular-buffer index of most recently stored pair
    active_old = .false.
    n_updates = 0
    n_rejects = 0
    n_resets  = 0
    conv_reason = -1

    ! project initial point into the feasible box
    do i = 1, n
       x(i) = min(max(x(i), l(i)), u(i))
    end do

    ! one gradient evaluation at the starting point
    call fg(n, x, f, g)

    ! initial projected-gradient norm for grtol test
    gnorm0 = proj_gnorm(n, x, g, l, u, BTOL)
    if (gnorm0 == 0.0d0) gnorm0 = 1.0d0

    info = 1
    iter = 0

    do k = 1, maxit
       iter = k

       ! active-set
       do i = 1, n
          active(i) = .false.
          if (x(i) <= l(i) + BTOL .and. g(i) > 0.0d0) active(i) = .true.
          if (x(i) >= u(i) - BTOL .and. g(i) < 0.0d0) active(i) = .true.
       end do

       ! reset L-BFGS history if the active set changed (mirrors TAO BQNLS)
       active_set_changed = any(active .neqv. active_old)
       if (active_set_changed .and. k > 1) then
          if (iprint > 0) write(*,'(A,I5)') &
               'BFGS: active set changed at iter ', k, ' -- resetting history'
          ncur   = 0
          newest = 0
          S = 0.0d0;  Y = 0.0d0;  rho = 0.0d0
          Bdiag = 1.0d0
          H0diag_stored = 1.0d0
          n_resets = n_resets + 1
       end if
       active_old = active

       ! projected-gradient norm on the free subspace
       gnorm = proj_gnorm(n, x, g, l, u, BTOL)

       if (iprint > 0) then
          write(*,'(A,I5,4(A,ES12.4))') &
               'BFGS iter =', k,                                        &
               '   f = ',           f,                                   &
               '   ||g|| = ',       gnorm,                               &
               '   ||g||/||g0|| =', gnorm / gnorm0,                      &
               '   ||g||/|f| = ',   gnorm / max(abs(f), 1.0d-300)
       end if

       ! --- three-way convergence test:
       !   gatol : ||g||        <= gatol
       !   grtol : ||g||/|f|   <= grtol
       !   gttol : ||g||/||g0|| <= gttol
       ! Convergence declared when ANY criterion is satisfied.
       if (gnorm <= gatol) then
          write(*,'(A,ES12.4,A,ES12.4)') &
               'BFGS: converged gatol ||g|| = ', gnorm, ' <= ', gatol
          info = 0;  conv_reason = 0;  gnorm_final = gnorm;  exit
       end if
       if (abs(f) > 0.0d0 .and. gnorm / abs(f) <= grtol) then
          write(*,'(A,ES12.4,A,ES12.4)') &
               'BFGS: converged grtol ||g||/|f| = ', gnorm/abs(f), &
               ' <= ', grtol
          info = 0;  conv_reason = 1;  gnorm_final = gnorm;  exit
       end if
       if (gnorm / gnorm0 <= gttol) then
          write(*,'(A,ES12.4,A,ES12.4)') &
               'BFGS: converged gttol ||g||/||g0|| = ', gnorm/gnorm0, &
               ' <= ', gttol
          info = 0;  conv_reason = 2;  gnorm_final = gnorm;  exit
       end if
       gnorm_final = gnorm   ! kept up to date each iteration in case
                              ! maxit or line-search-failure exit occurs

       ! --- L-BFGS two-loop recursion (Nocedal & Wright, Alg. 7.4)
       !
       ! All dot products and updates are masked to the FREE subspace.
       ! Active components of q/d are re-zeroed after every update step
       ! to prevent stale S/Y components at those indices from leaking
       ! back in.  rho values are also computed on the free subspace
       ! (see storage below).

       q = g
       do i = 1, n
          if (active(i)) q(i) = 0.0d0
       end do

       ! backward pass
       do i = 1, ncur
          ! index into circular buffer: newest pair first
          col = mod(newest - i + m, m) + 1
          alpha_ls(i) = rho(col) * dot_product(S(:,col), q)
          q = q - alpha_ls(i) * Y(:,col)
          do j = 1, n          ! re-mask: Y may be non-zero at active indices
             if (active(j)) q(j) = 0.0d0
          end do
       end do

       ! --- H0 initialization: SCALAR (h0_type=0) or DIAGONAL (h0_type=1)
       if (h0_type == 0) then
          ! SCALAR: classic Nocedal gamma = (s^T y)/(y^T y) using the
          ! most recent accepted pair, restricted to the free subspace.
          if (ncur > 0) then
             col = newest
             sty_free = 0.0d0;  yty_free = 0.0d0
             do i = 1, n
                if (.not. active(i)) then
                   sty_free = sty_free + S(i,col) * Y(i,col)
                   yty_free = yty_free + Y(i,col) * Y(i,col)
                end if
             end do
             H0scale = sty_free / max(yty_free, 1.0d-300)
             H0scale = max(H0scale, 1.0d-10)
          else
             H0scale = 1.0d0
          end if
          do i = 1, n
             if (active(i)) then
                H0diag(i) = 1.0d0
             else
                H0diag(i) = H0scale
             end if
          end do
       else
          ! DIAGONAL: per-component Broyden-class scaling. H0diag_stored
          ! is the fully-computed (theta-blend + sigma-rescale) inverse
          ! diagonal, updated once per accepted step below -- just use
          ! it directly here.
          do i = 1, n
             if (active(i)) then
                H0diag(i) = 1.0d0
             else
                H0diag(i) = H0diag_stored(i)
             end if
          end do
       end if
       d = H0diag * q

       ! forward pass
       do i = ncur, 1, -1
          col = mod(newest - i + m, m) + 1
          beta_ls = rho(col) * dot_product(Y(:,col), d)
          d = d + S(:,col) * (alpha_ls(i) - beta_ls)
          do j = 1, n          ! re-mask after S update
             if (active(j)) d(j) = 0.0d0
          end do
       end do

       d = -d   ! search direction

       ! safeguard: fall back to scaled steepest descent on the free set
       if (dot_product(d, g) > -1.0d-12) then
          d = 0.0d0
          do i = 1, n
             if (.not. active(i)) d(i) = -H0diag(i) * g(i)
          end do
       end if

       ! --- backtracking Armijo line search (cheap: b2mn_step only)
       xold = x
       fold = f
       gold = g

       ! Initial trial step length. On the VERY FIRST outer iteration
       ! (no curvature pairs exist yet, so d = -g exactly regardless of
       ! h0_type), use a special heuristic rather than the plain 
       ! configured stepinit: the classic Nocedal & Wright initial-step formula
       !     alpha_0 = 2*|f_0| / ||g_0||^2
       ! (implicitly assuming an unknown-but-optimistic target f_low=0
       ! for the expected decrease). Without this, a badly-scaled problem's 
       ! first step is a full, unscaled H0=I Newton-like step that can be 
       ! enormous and immediately drive every variable to its bound.
       ! For k>1 (once curvature/H0 scaling is available), the usual
       ! configured alpha_init is used as before.
       if (k == 1 .and. ncur == 0) then
          alpha = min(alpha_init, 2.0d0*abs(f) / max(gnorm**2, 1.0d-300))
          if (iprint > 0) write(*,'(A,ES12.4)') &
               'BFGS: first-iteration initial step length = ', alpha
       else
          alpha = alpha_init   ! matches TAO's -tao_ls_stepinit (default 1.0)
       end if

       do
          do i = 1, n
             x(i) = min(max(xold(i) + alpha*d(i), l(i)), u(i))
          end do
          call f_only(n, x, f)

          pred = dot_product(gold, x - xold)
          if (f <= fold + C1 * pred) exit

          alpha = alpha * RHOLS
          if (alpha < ALPHA_MIN) then
             info = 2;  exit
          end if
       end do

       if (info == 2) exit

       ! --- expensive gradient evaluation at the accepted point
       call fg(n, x, f, g)

       ! --- store new (s, y) pair; rho computed on free subspace only
       ! (variables active at the *new* point are excluded so stale
       ! curvature from the bound doesn't enter future directions)
       sty_free = 0.0d0
       do i = 1, n
          if (.not. active(i)) &
             sty_free = sty_free + (x(i)-xold(i)) * (g(i)-gold(i))
       end do

       if (sty_free > CURV_TOL) then
          newest     = mod(newest, m) + 1
          ncur       = min(ncur + 1, m)
          S(:,newest) = x - xold
          Y(:,newest) = g - gold
          rho(newest) = 1.0d0 / sty_free
          n_updates   = n_updates + 1

          ! --- diagonal Broyden update: ONLY needed/computed when
          ! h0_type=1 (DIAGONAL). Skipped entirely for h0_type=0
          ! (SCALAR) so the scalar path has zero extra computation and
          ! is unaffected by any of the diagonal-specific machinery.
          ! DIAGONAL
          ! This is a convex combination of the diagonal-restricted
          ! BFGS and DFP updates:
          !   Bdiag_new = (1-theta)*B_BFGS + theta*B_DFP
          ! with theta=0 reducing to pure BFGS, which has TWO stages:
          !
          ! STAGE 1 -- per-component theta-blended Broyden update of
          ! the DIRECT Hessian diagonal Bdiag:
          !   sBs   = s^T Bdiag_old s          (restricted to free set)
          !   stDs  = max(sBs, RESCALE_TOL)    (PETSc's ldb->tol, 1e-8)
          !   bfgs_i = -(Bdiag_old(i)*s(i))^2 / stDs
          !   dfp_i  = stDs/sty * y(i)^2 - 2*y(i)*Bdiag_old(i)*s(i)
          !   Bdiag_new(i) = Bdiag_old(i) + y(i)^2/sty
          !                  + (1-theta)*bfgs_i + theta*dfp_i/sty
          ! (the "+y^2/sty" term is added in FULL regardless of theta)
          !
          ! STAGE 2 -- reciprocal + abs() + a GLOBAL scalar rescale
          ! "sigma" applied on top:
          !   H0(i) = abs(1/Bdiag_new(i))
          !   sigma = sty / sum_i(H0(i)*y(i)^2)     [PETSc defaults:
          !           alpha=1, beta=0.5, sigma_hist=1 -- i.e. only the
          !           MOST RECENT pair is used; general sigma_hist>1
          !           history windowing is not implemented here]
          !   H0(i) = H0(i) * sigma   (only if sigma is finite & > 0)
          !
          ! STAGE 3 -- convex blend with the previous H0 via rescale_rho:
          !   H0_final(i) = (1-rescale_rho)*H0_old(i) + rescale_rho*H0(i)
          !
          ! BDIAG_MIN/MAX below are defensive clamps, kept as a generous
          ! safety net.
          if (h0_type == 1) then
             sBs = 0.0d0
             do i = 1, n
                if (.not. active(i)) &
                   sBs = sBs + Bdiag(i) * (S(i,newest))**2
             end do
             sBs = max(sBs, RESCALE_TOL)

             do i = 1, n
                if (.not. active(i)) then
                   bfgs_i = -(Bdiag(i)*S(i,newest))**2 / sBs
                   dfp_i  = sBs/sty_free * (Y(i,newest))**2 &
                          - 2.0d0*Y(i,newest)*Bdiag(i)*S(i,newest)
                   Bdiag(i) = Bdiag(i) + (Y(i,newest))**2/sty_free &
                            + (1.0d0 - theta_c)*bfgs_i &
                            + theta_c*dfp_i/sty_free
                   Bdiag(i) = min(max(Bdiag(i), BDIAG_MIN), BDIAG_MAX)
                end if
             end do

             ! --- stage 2: reciprocal, abs
             do i = 1, n
                if (.not. active(i)) then
                   H0new(i) = abs(1.0d0 / Bdiag(i))
                end if
             end do

             ! --- stage 2b: global sigma rescale (sigma_hist=1 default:
             ! uses only the pair just stored)
             yHy = 0.0d0
             do i = 1, n
                if (.not. active(i)) yHy = yHy + H0new(i) * (Y(i,newest))**2
             end do
             if (yHy > 1.0d-300) then
                sigma_scale = sty_free / yHy
                if (sigma_scale == sigma_scale .and. &   ! reject NaN
                    abs(sigma_scale) < huge(1.0d0) .and. &  ! reject Inf
                    sigma_scale > 0.0d0) then
                   do i = 1, n
                      if (.not. active(i)) H0new(i) = H0new(i) * sigma_scale
                   end do
                end if
             end if

             ! --- stage 3: rescale_rho blend (default 1 = full replace)
             do i = 1, n
                if (.not. active(i)) then
                   if (rescale_rho_c >= 1.0d0) then
                      H0diag_stored(i) = H0new(i)
                   else
                      H0diag_stored(i) = (1.0d0 - rescale_rho_c)*H0diag_stored(i) &
                                       + rescale_rho_c*H0new(i)
                   end if
                end if
             end do
          end if   ! h0_type == 1
       else
          n_rejects = n_rejects + 1
       end if


    end do

    deallocate(S, Y, rho, g, gold, xold, d, q, alpha_ls, active, active_old, &
               Bdiag, H0diag, H0new, H0diag_stored)
  end subroutine bfgs_bound_solve


  !---------------------------------------------------------------------
  ! proj_gnorm: norm of g restricted to the free subspace
  !---------------------------------------------------------------------
  pure function proj_gnorm(n, x, g, l, u, btol) result(nrm)
    integer,          intent(in) :: n
    double precision, intent(in) :: x(n), g(n), l(n), u(n), btol
    double precision :: nrm
    integer :: i
    nrm = 0.0d0
    do i = 1, n
       if (x(i) <= l(i) + btol .and. g(i) > 0.0d0) cycle
       if (x(i) >= u(i) - btol .and. g(i) < 0.0d0) cycle
       nrm = nrm + g(i)*g(i)
    end do
    nrm = sqrt(nrm)
  end function proj_gnorm

end module bfgs_bound_mod


!=======================================================================
! Main optimisation program
!=======================================================================
program b2optim_bfgs
  use bfgs_bound_mod ! IGNORE
  use b2mod_par_opt_diff
  use b2mod_main_diff  &
     , only : b2mn_init_diff, b2mn_step, b2mn_step_diff, b2mn_fin_diff
  use b2mod_ad_diff    &
     , only : nncf
  use b2us_data_diff
  implicit none

  ! plain arrays replace PETSc Vec objects
  real(kind=R8), allocatable :: X(:)    ! current solution (rescaled)
  real(kind=R8), allocatable :: X_L(:)  ! lower bounds    (rescaled)
  real(kind=R8), allocatable :: X_U(:)  ! upper bounds    (rescaled)

  ! sentinel value for "no bound" (replaces PETSC_INFINITY)
  real(kind=R8), parameter :: INF_BOUND = 1.0e30_R8

  integer, save :: iter  = 0
  integer, save :: filen = 0
  integer, save :: ntim  = 1

  integer :: bfgs_iter, bfgs_info, bfgs_conv_reason, ipar
  integer :: bfgs_n_updates, bfgs_n_rejects, bfgs_n_resets
  real(kind=R8) :: bfgs_gnorm_final
  logical :: streql
  external streql

  ! ---- initialise B2 physics ----------------------------------------
  flag_optim = .true.
  call b2mn_init_diff()
  call ipgeti('b2mndr_ntim', ntim)
  par_opt_phys = 0.0_R8

#ifdef TGT
  call xertst(npar_opt .le. nbdirsmax, &
              'Increase size of nbdirsmax in b2mod_diffsizes.F')
  call set_tgt_perturbation(switchd)
#endif

#ifdef ADJ
  call print_adj_parameters()
  par_opt_physb = 0.0_R8
#endif

  ! ---- build plain-array problem ------------------------------------
  call InitializeProblem(npar_opt)

  ! ---- call the BFGS solver -----------------------------------------
  ! tol_opt and maxiter come from b2mod_par_opt_diff (same as TAO path)
  call xertst(tol_opt  > 0.0_R8, 'faulty internal parameter tol_opt')
  call xertst(maxiter  > 0,       'faulty internal parameter maxiter')

  ! ---- L-BFGS memory size (mirrors TAO's -tao_bqnls_mat_lmvm_hist_size)
  ! Default = 1 to match TAO's out-of-the-box TAOBQNLS behaviour.
  ! Increase (e.g. 5, 10) for potentially faster convergence at the
  ! cost of more memory; useful when npar_opt is small.
  call xertst(0.lt.lbfgs_memsize, 'wrong parameter lbfgs_memsize')
  write(*,'(A,I0)') ' b2optim_bfgs: L-BFGS memory size = ', lbfgs_memsize

  ! ---- initial-Hessian (J0) scaling scheme
  ! 0 = SCALAR (classic Nocedal gamma from the most recent pair)
  ! 1 = DIAGONAL (per-component Broyden-class scaling; see theta below)
  !
  ! *** On at least one real test case, SCALAR (0) converged in 18
  ! iterations, matching TAO's 16 closely, while DIAGONAL at theta=0
  ! took 22 -- these are genuinely different algorithms and which one
  ! is faster is problem-dependent, not a strict "one is more correct"
  ! situation. Default here is SCALAR (0) since it was the best match
  ! found so far; try DIAGONAL (1) if you want to experiment, e.g. for
  ! problems with many parameters of very different scales.
  call xertst((lbfgs_h0_type.eq.0) .or. (lbfgs_h0_type.eq.1), &
              'wrong parameter lbfgs_h0_type')
  write(*,'(A,I0,A)') ' b2optim_bfgs: H0 scaling type = ', lbfgs_h0_type, &
       merge(' (SCALAR)  ', ' (DIAGONAL)', lbfgs_h0_type == 0)

  ! ---- L-BFGS diagonal J0 convex mixing factor
  ! Mirrors TAO's -tao_bqnls_mat_lmvm_theta (TAO's own default is 0.125).
  !
  ! *** IMPORTANT CAVEAT ***
  ! Testing found that our simplified diagonal-DFP term can enter a
  ! positive-feedback loop under repeated recursive reuse (each update
  ! feeds its own output back into the next one via sBs), causing the
  ! solver to stagnate for theta > 0 on at least one test problem, even
  ! with a growth-ratio limiter applied. PETSc's actual implementation
  ! has additional stabilization (alpha/beta/rho/mu/nu parameters) that
  ! are not publicly documented in enough detail to replicate exactly,
  ! so this theta blend should be treated as EXPERIMENTAL.
  !
  ! Default is therefore kept at 0 (pure BFGS diagonal scaling), which
  ! was validated to work well and even outperform TAO's 16-iteration
  ! result in some configurations. Only set this above 0 if you want to
  ! experiment, and verify convergence carefully (check BFGS.OUT status
  ! and updates/rejects/resets) before trusting the result.
  call xertst((0.0_R8.le.lbfgs_theta).and.(lbfgs_theta.le.1.0_R8), &
              'wrong parameter lbfgs_theta')
  write(*,'(A,F8.4)') ' b2optim_bfgs: L-BFGS diagonal theta = ', lbfgs_theta

  ! ---- rescale_rho: convex blend factor between the old and newly
  ! rescaled H0 diagonal (PETSc's -mat_lmvm_rho, default 1.0 = fully
  ! replace with the new value each update, i.e. a no-op blend). Only
  ! meaningful when h0_type=1 (DIAGONAL).
  call xertst((0.0_R8.le.lbfgs_rescale_rho).and.(lbfgs_rescale_rho.le.1.0_R8), &
              'wrong parameter lbfgs_rescale_rho')
  write(*,'(A,F8.4)') ' b2optim_bfgs: rescale rho = ', lbfgs_rescale_rho

  ! ---- initial line-search step length
  ! Mirrors TAO's -tao_ls_stepinit (TAO's own default is 1.0). This is
  ! the trial step length alpha tried FIRST at every outer iteration,
  ! before any backtracking. Lowering it can help if the very first
  ! full Newton-like step tends to overshoot badly (e.g. very
  ! nonlinear/stiff physics); raising it above 1 is unusual but
  ! allowed.
  call xertst((0.0_R8.lt.lbfgs_stepinit).and.(lbfgs_stepinit.le.1.0_R8), &
              'wrong parameter lbfgs_stepinit')
  write(*,'(A,F8.4)') ' b2optim_bfgs: line search initial step = ', lbfgs_stepinit

  ! ---- call the L-BFGS solver ---------------------------------------
  ! tol_opt for all three tolerances matches:
  !   call TaoSetTolerances(tao, tol_opt, tol_opt, tol_opt, ierr)
  call bfgs_bound_solve(npar_opt, X, X_L, X_U,               &
                        FormFunction,                        &  ! cheap  (b2mn_step)
                        FormFunctionGradient,                &  ! expensive (b2mn_step_diff)
                        maxiter,                             &
                        real(tol_opt,8),                     &  ! gatol
                        real(tol_opt,8),                     &  ! grtol
                        real(tol_opt,8),                     &  ! gttol
                        lbfgs_memsize,                       &  ! L-BFGS history
                        lbfgs_h0_type,                       &  ! H0 scaling type
                        real(lbfgs_theta,8),                 &  ! diagonal mixing
                        real(lbfgs_rescale_rho,8),           &  ! rescale rho blend
                        real(lbfgs_stepinit,8),              &  ! line search step init
                        1,                                   &  ! iprint
                        bfgs_iter, bfgs_info, bfgs_conv_reason, &
                        bfgs_gnorm_final,                    &
                        bfgs_n_updates, bfgs_n_rejects, bfgs_n_resets)

  ! ---- report -------------------------------------------------------
  call WriteResults(bfgs_iter, bfgs_info, bfgs_conv_reason,  &
                     bfgs_gnorm_final,                       &
                     bfgs_n_updates, bfgs_n_rejects, bfgs_n_resets)

  ! ---- finalise -----------------------------------------------------
  call DestroyProblem()
  call b2mn_fin_diff()
  deallocate(par_opt_phys)
  deallocate(xold)
  deallocate(xnew)
  deallocate(xmult)
  deallocate(par_opt_physdiff)
  stop 'b2optim_bfgs'

contains

  !=====================================================================
  ! InitializeProblem
  ! Allocate and fill X, X_L, X_U from b2mod_par_opt_diff arrays,
  ! mirroring what VecSetValue / VecAssemblyEnd did in b2optim_tao.
  !=====================================================================
  subroutine InitializeProblem(npar)
    use b2mod_par_opt_diff, only : x0, xl, xu, par_rescale
    integer, intent(in) :: npar

    allocate(X(npar), X_L(npar), X_U(npar))

    do ipar = 1, npar

       X(ipar)   = x0(ipar) / par_rescale(ipar)
       xold(ipar) = x0(ipar) / par_rescale(ipar)

       if (xl(ipar) < -inf_opt) then
          write(*,*) 'BFGS: warning, X_L(', ipar, ') set to -INF_BOUND'
          X_L(ipar) = -INF_BOUND
       else
          X_L(ipar) = xl(ipar) / par_rescale(ipar)
       end if

       if (xu(ipar) > inf_opt) then
          write(*,*) 'BFGS: warning, X_U(', ipar, ') set to +INF_BOUND'
          X_U(ipar) = INF_BOUND
       else
          X_U(ipar) = xu(ipar) / par_rescale(ipar)
       end if

    end do

    write(*,*) 'BFGS SET X0'
    write(*,'(*(ES14.6,2X))') X(1:npar)
    write(*,*) 'BFGS SET X_L'
    write(*,'(*(ES14.6,2X))') X_L(1:npar)
    write(*,*) 'BFGS SET X_U'
    write(*,'(*(ES14.6,2X))') X_U(1:npar)

  end subroutine InitializeProblem


  !=====================================================================
  ! DestroyProblem
  !=====================================================================
  subroutine DestroyProblem()
    if (allocated(X))   deallocate(X)
    if (allocated(X_L)) deallocate(X_L)
    if (allocated(X_U)) deallocate(X_U)
  end subroutine DestroyProblem


  !=====================================================================
  ! WriteResults
  ! Replaces TaoView / PetscViewerASCIIOpen output.
  !=====================================================================
  subroutine WriteResults(it, info, conv_reason, gnorm_final, &
                           n_updates, n_rejects, n_resets)
    integer,          intent(in) :: it, info, conv_reason
    double precision, intent(in) :: gnorm_final
    integer,          intent(in) :: n_updates, n_rejects, n_resets
    integer :: iu
    character(len=80) :: status_line, reason_line

    iu = 99
    open(unit=iu, file='BFGS.OUT', status='replace')
    write(iu,'(A)')    '========================================='
    write(iu,'(A)')    ' b2optim_bfgs  --  BFGS optimisation'
    write(iu,'(A)')    '========================================='
    write(iu,'(A,I0)') ' Total iterations : ', it
    write(iu,'(A,ES12.4,A,ES12.4,A,ES12.4)') &
         ' Tolerances (gatol/grtol/gttol): ', &
         tol_opt, ' / ', tol_opt, ' / ', tol_opt
    write(iu,'(A,ES14.6)') ' Final projected gradient norm ||g|| : ', gnorm_final

    ! --- explicit convergence-criterion report ---------------------
    ! Reports EXACTLY which of the three stopping conditions was met,
    ! or which failure mode caused termination -- not just "converged".
    select case (info)
    case (0)
      status_line = 'CONVERGED'
      select case (conv_reason)
      case (0)
        reason_line = 'gatol criterion met:  ||g|| <= gatol'
      case (1)
        reason_line = 'grtol criterion met:  ||g||/|f| <= grtol'
      case (2)
        reason_line = 'gttol criterion met:  ||g||/||g0|| <= gttol'
      case default
        reason_line = '(unknown -- should not happen)'
      end select
    case (1)
      status_line  = 'STOPPED'
      reason_line  = 'maximum iteration count reached before any ' // &
                     'tolerance was satisfied'
    case (2)
      status_line  = 'FAILED'
      reason_line  = 'line search step length collapsed below the ' // &
                     'minimum allowed step (no acceptable decrease found)'
    case default
      status_line  = 'UNKNOWN'
      reason_line  = '(unrecognized info code -- should not happen)'
    end select

    write(iu,'(A)') ' '
    write(iu,'(A,A)') ' Status             : ', trim(status_line)
    write(iu,'(A,A)') ' Convergence reason : ', trim(reason_line)

    write(iu,'(A)') ' '
    write(iu,'(A)') ' Final solution (rescaled):'
    write(iu,'(*(ES14.6,2X))') X(1:npar_opt)
    write(iu,'(A)') ' Final solution (physical):'
    write(iu,'(*(ES14.6,2X))') X(1:npar_opt) * par_rescale(1:npar_opt)

    ! --- L-BFGS matrix summary, mirroring TAO's "Mat Object" view:
    !   Mat Object: (tao_bqnls_) ... type: lmvmbfgs
    !     Max. storage / Used storage / Number of updates / rejects / resets
    write(iu,'(A)') ' '
    write(iu,'(A)') ' L-BFGS matrix summary:'
    write(iu,'(A,I0)') '   Max. storage        : ', lbfgs_memsize
    write(iu,'(A,I0)') '   Used storage        : ', min(n_updates, lbfgs_memsize)
    write(iu,'(A,I0)') '   Number of updates   : ', n_updates
    write(iu,'(A,I0)') '   Number of rejects   : ', n_rejects
    write(iu,'(A,I0)') '   Number of resets    : ', n_resets
    write(iu,'(A,F8.4)') '   Rescale theta       : ', lbfgs_theta
    write(iu,'(A,F8.4)') '   Line search stepinit: ', lbfgs_stepinit
    close(iu)

    ! also echo to stdout
    write(*,*) '========================================='
    write(*,*) ' b2optim_bfgs  --  BFGS optimisation'
    write(*,*) ' Total iterations : ', it
    write(*,*) ' Status             : ', trim(status_line)
    write(*,*) ' Convergence reason : ', trim(reason_line)
    write(*,'(A,I0,A,I0,A,I0)') &
         ' L-BFGS updates/rejects/resets : ', n_updates, ' / ', &
         n_rejects, ' / ', n_resets
  end subroutine WriteResults


  !=====================================================================
  ! FormFunction
  !
  ! Cheap (forward-only) evaluator called by the BFGS line search.
  ! Wraps b2mn_step -- no adjoint/tangent solve.
  ! Signature matches func_if abstract interface.
  !=====================================================================
  subroutine FormFunction(n, x_v, F)
    use b2mod_par_opt_diff, only : sigma, mean, par_rescale
    use b2mod_ad_diff,      only : primal_iterations, primal_res
    integer,          intent(in)  :: n
    double precision, intent(in)  :: x_v(n)
    double precision, intent(out) :: F

    real(kind=R8) :: j(nncf)
    integer :: ipar, isigma, imean, ishift, icorr
    character*3 :: str

    ! pack x_v into the B2 parameter arrays (same logic as FormFunction
    ! in b2optim_tao)
    call reset_drifts_params(n, x_v)

    do ipar = 1, npar_opt - nsigma_opt - nmean_opt - nshift_opt - ncorr_opt
       par_opt_phys(ipar) = xold(ipar) * par_rescale(ipar)
       xnew(ipar)         = x_v(ipar)  * par_rescale(ipar)
       write(str,"(I1)") ipar
       if (ipar >= 10) write(str,"(I2)") ipar
       write(*,*) 'BFGS: eval_F with x',trim(str),'= ', xnew(ipar)
    end do

    if (cftype(1) == 6) then
       isigma = npar_opt - nsigma_opt - nmean_opt - nshift_opt - ncorr_opt + 1
       do ipar = 1, nsigma
          if (sigma_opt(ipar)) then
             sigma(ipar) = x_v(isigma) * par_rescale(isigma)
             write(str,"(I1)") isigma
             if (isigma >= 10) write(str,"(I2)") isigma
             write(*,*) 'BFGS: eval_F with x',trim(str),'= ', sigma(ipar)
             isigma = isigma + 1
          end if
       end do
       imean = npar_opt - nmean_opt - nshift_opt - ncorr_opt + 1
       do ipar = 1, nmean
          if (mean_opt(ipar)) then
             mean(ipar) = x_v(imean) * par_rescale(imean)
             write(str,"(I1)") imean
             if (imean >= 10) write(str,"(I2)") imean
             write(*,*) 'BFGS: eval_F with x',trim(str),'= ', mean(ipar)
             imean = imean + 1
          end if
       end do
    end if

    ishift = npar_opt - nshift_opt - ncorr_opt + 1
    do ipar = 1, nshift
       if (shiftopt(ipar)) then
          shift(ipar) = x_v(ishift) * par_rescale(ishift)
          write(str,"(I1)") ishift
          if (ishift >= 10) write(str,"(I2)") ishift
          write(*,*) 'BFGS: eval_F with x',trim(str),'= ', shift(ipar)
          ishift = ishift + 1
       end if
    end do

    if (cftype(1) == 6) then
       icorr = npar_opt - ncorr_opt + 1
       do ipar = 1, ncf
          if (corr_opt(ipar)) then
             corr_length(ipar) = x_v(icorr) * par_rescale(icorr)
             write(str,"(I1)") icorr
             if (icorr >= 10) write(str,"(I2)") icorr
             write(*,*) 'BFGS: eval_F with x',trim(str),'= ', corr_length(ipar)
             icorr = icorr + 1
          end if
       end do
    end if

    ! forward-only solve (cheap)
    call b2mn_step(j)
    xold(1:npar_opt) = x_v(1:npar_opt)
    F = real(j(1), 8)

    write(*,*) 'BFGS COST FUNCTION:', F
    write(*,*) 'BFGS PRIMAL ITERATIONS:', primal_iterations
    write(*,*) 'BFGS PRIMAL RESIDUAL:',   primal_res

  end subroutine FormFunction


  !=====================================================================
  ! FormFunctionGradient
  !
  ! Expensive evaluator (forward + adjoint/tangent) called once per
  ! accepted BFGS iterate.  Wraps b2mn_step_diff.
  ! Signature matches func_grad_if abstract interface.
  !=====================================================================
  subroutine FormFunctionGradient(n, x_v, F, g_v)
    use b2us_io_diff,       only : write_b2fstate
    use b2mod_version,       only : newversion, cfverw
    use b2mod_b2cmpa
    use b2mod_par_opt_diff, only : par_rescale, sigma, mean
    use b2mod_ad_diff,      only : primal_iterations, gradient_iterations, &
                                    primal_res, gradient_res
    integer,          intent(in)  :: n
    double precision, intent(in)  :: x_v(n)
    double precision, intent(out) :: F
    double precision, intent(out) :: g_v(n)

    real(kind=R8) :: j(nncf), jdiff(nncf), gradd(npar_opt)
    integer :: ipar, isigma, imean, ishift, icorr, idum(0:2)
    character*3  :: str
    character(22) :: opt_state_name
    character*120 :: label
    logical :: write_state

    call reset_drifts_params(n, x_v)

    do ipar = 1, npar_opt - nsigma_opt - nmean_opt - nshift_opt - ncorr_opt
       par_opt_phys(ipar) = xold(ipar) * par_rescale(ipar)
       xnew(ipar)         = x_v(ipar)  * par_rescale(ipar)
       write(str,"(I1)") ipar
       if (ipar >= 10) write(str,"(I2)") ipar
       write(*,*) 'BFGS: eval_F_grad_F with x',trim(str),'= ', xnew(ipar)
    end do

    if (cftype(1) == 6) then
       isigma = npar_opt - nsigma_opt - nmean_opt - nshift_opt - ncorr_opt + 1
       do ipar = 1, nsigma
          if (sigma_opt(ipar)) then
             sigma(ipar) = x_v(isigma) * par_rescale(isigma)
             write(str,"(I1)") isigma
             if (isigma >= 10) write(str,"(I2)") isigma
             write(*,*) 'BFGS: eval_F_grad_F with x',trim(str),'= ', sigma(ipar)
             isigma = isigma + 1
          end if
       end do
       imean = npar_opt - nmean_opt - nshift_opt - ncorr_opt + 1
       do ipar = 1, nmean
          if (mean_opt(ipar)) then
             mean(ipar) = x_v(imean) * par_rescale(imean)
             write(str,"(I1)") imean
             if (imean >= 10) write(str,"(I2)") imean
             write(*,*) 'BFGS: eval_F_grad_F with x',trim(str),'= ', mean(ipar)
             imean = imean + 1
          end if
       end do
    end if

    ishift = npar_opt - nshift_opt - ncorr_opt + 1
    do ipar = 1, nshift
       if (shiftopt(ipar)) then
          shift(ipar) = x_v(ishift) * par_rescale(ishift)
          write(str,"(I1)") ishift
          if (ishift >= 10) write(str,"(I2)") ishift
          write(*,*) 'BFGS: eval_F_grad_F with x',trim(str),'= ', shift(ipar)
          ishift = ishift + 1
       end if
    end do

    if (cftype(1) == 6) then
       icorr = npar_opt - ncorr_opt + 1
       do ipar = 1, ncf
          if (corr_opt(ipar)) then
             corr_length(ipar) = x_v(icorr) * par_rescale(icorr)
             write(str,"(I1)") icorr
             if (icorr >= 10) write(str,"(I2)") icorr
             write(*,*) 'BFGS: eval_F_grad_F with x',trim(str),'= ', corr_length(ipar)
             icorr = icorr + 1
          end if
       end do
    end if

    ! forward + gradient solve
    ! neff: exclude sigma/mean/shift/corr parameters from the tangent/adjoint
    ! run (same reasoning as in FormFunctionGradient in b2optim_tao)
    call b2mn_step_diff(j, jdiff)
    xold(1:npar_opt) = x_v(1:npar_opt)
    F = real(j(1), 8)

#ifdef TGT
    do ipar = 1, npar_opt
       g_v(ipar) = real(jdiff(1)*par_rescale(ipar), 8)
#endif
#ifdef ADJ
    call set_adj_gradient(npar_opt, gradd, switchdiff)
    do ipar = 1, npar_opt
       g_v(ipar) = real(gradd(ipar)*par_rescale(ipar), 8)
#endif
       write(str,"(I1)") ipar
       if (ipar >= 10) write(str,"(I2)") ipar
       write(*,*) 'BFGS GRAD of x',trim(str),':', g_v(ipar)
    end do

    write(*,*) 'BFGS COST FUNCTION:',        F
    write(*,*) 'BFGS GRADIENT NORM:',         norm2(g_v(1:npar_opt))
    write(*,*) 'BFGS PRIMAL ITERATIONS:',     primal_iterations
    write(*,*) 'BFGS GRADIENT ITERATIONS:',   gradient_iterations
    write(*,*) 'BFGS PRIMAL RESIDUAL:',       primal_res
    write(*,*) 'BFGS GRADIENT RESIDUAL:',     gradient_res

    ! --- optional intermediate state dump (mirrors b2optim_tao logic)
    write_state = .false.
    if (switch%b2optim_save_states > 0) then
       write_state = mod(iter, switch%b2optim_save_states) == 0
    end if

    if (write_state) then
       write(*,*) 'Saving intermediate optimization state'
       write(opt_state_name,'(a14,i4.4)') 'b2fstate_optim.', filen
       call cfopen(99, trim(opt_state_name), 'new', 'un*formatted')
       idum(0) = mpg%nCv
       idum(1) = mpg%nFc
       idum(2) = state%ns
       call cfverw(99, newversion)
       call cfwuin(99, 3, idum, 'nCv,nFc,ns')
       write(label,'(a46,i4)') 'b2optim_bfgs intermediate optimization state ', filen
       call cfwuch(99, 120, label, 'label')
       call b2wuzd(99, newversion, state%ns, zamin, zamax, zn, am)
       call write_b2fstate(99, mpg%nCv, mpg%nFc, state%ns, state)
       close(99)
       filen = filen + 1
    end if

    iter = iter + 1

  end subroutine FormFunctionGradient


  !=====================================================================
  ! reset_drifts_params
  !
  ! Identical in purpose and logic to the same-named subroutine in
  ! b2optim_tao.  Operates on plain arrays instead of PETSc Vec.
  !=====================================================================
  subroutine reset_drifts_params(n, x_v)
    use b2mod_facdrift_exb_diff, only : facdrift_scalar, fac_exb_scalar
    integer,          intent(in) :: n
    double precision, intent(in) :: x_v(n)

    logical :: sameX
    integer :: ipar
    real(R8) :: nn

    ! check whether x_v is the same as the last accepted point
    sameX = .true.
    do ipar = 1, npar_opt
       if (abs(xold(ipar) - x_v(ipar)) > 1.0e-4_R8) sameX = .false.
    end do

    if (.not. sameX .and. &
        (switch%facExB_start > 0.0_R8 .or. switch%facdrift_start > 0.0_R8)) then
       ! reset drift fractions and rebuild the ramp increment
       fac_exb_scalar  = min(switch%b2optim_reset_drift, switch%facExB_start)
       facdrift_scalar = min(switch%b2optim_reset_drift, switch%facdrift_start)
       nn = real(ntim, R8) / switch%b2optim_reset_drift_iter
       switch%facExB_inc   = (1.0_R8 / fac_exb_scalar )**(1.0_R8 / nn)
       switch%facdrift_inc = (1.0_R8 / facdrift_scalar)**(1.0_R8 / nn)
       write(*,*) ' b2optim_bfgs: resetting drifts to',  switch%b2optim_reset_drift
       write(*,*) ' b2optim_bfgs: drift_inc set to',     switch%facExB_inc
       write(*,*) ' b2optim_bfgs: drifts back at 100% in ', nint(nn), ' iterations'
    end if

    if (.not. sameX .and. switch%b2optim_reset_param_iter > 1) then
       recalc_params = .true.
       do ipar = 1, npar_opt
          nn = real(ntim, R8) / switch%b2optim_reset_param_iter
          xmult(ipar) = (x_v(ipar) / xold(ipar))**(1.0_R8 / nn)
       end do
    else
       recalc_params = .false.
       xold(1:npar_opt) = x_v(1:npar_opt)
    end if

  end subroutine reset_drifts_params

end program b2optim_bfgs
