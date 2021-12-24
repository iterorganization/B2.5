      module taomodule
#include <petsc/finclude/petsctao.h>
      use petsctao ! IGNORE

      Vec X,X_L,X_U
      Mat Hess
      end module taomodule

      program b2optim_petsc
      use taomodule ! IGNORE
      use b2mod_par_opt_diff
      use b2mod_main_diff &
      , only : b2mn_init_d, b2mn_step, b2mn_step_d, b2mn_fin_d
      use b2mod_ad_diff &
      , only : nncf
      use b2mod_switches_diff
      implicit none

      PetscErrorCode       ierr
      Tao                  tao
      KSP                  ksp
      PC                   pc
      PetscReal            tol
      PetscViewer          viewer

      integer :: ncon, nele_jac, ipar
      logical :: streql, hessian
      external streql
      type(switches_diffv), save :: switchdiff

      call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
      if (ierr .ne. 0) then
         print*,'Unable to initialize PETSc'
         stop
      endif

!      call PetscPrintf(PETSC_COMM_SELF,'Solution should be f(1,1)=-2\n',ierr);CHKERRA(ierr)

      ! Allocate and initialize par_opt variables to be used in B2.5
      flag_optim  = .true.
      call b2mn_init_d
      par_opt_phys = 0.0_R8
!     Initialize derivatives of estimated parameters
#ifdef TGT
      call xertst(npar_opt.le.nbdirsmax, 'Increase size of nbdirsmax in diffsizes.F')
      call set_tgt_perturbation(switchdiff)
#endif

#ifdef ADJ
      par_opt_physb = 0.0_R8
#endif

!     For the moment only steepest descent is available when specifiyng the Hessian
      hessian = .false.
      if (streql('exact',hessian_approximation)) then
        hessian = .true.
        write (*,*) 'PETSC-TAO, WARNING: setting Hessian as identity matrix for Steepest Descent method'
      endif

      call InitializeProblem(npar_opt,ierr);CHKERRA(ierr)

      call TaoCreate(PETSC_COMM_SELF,tao,ierr);CHKERRA(ierr)
!     Here better coding can be done to have different methods?
      if (hessian) then
        call TaoSetType(tao,TAOTRON,ierr);CHKERRA(ierr) 
      else
        call TaoSetType(tao,TAOBQNLS,ierr);CHKERRA(ierr) 
      endif
      call TaoSetInitialVector(tao,X,ierr);CHKERRA(ierr)
      call TaoSetVariableBounds(tao,X_L,X_U,ierr);CHKERRA(ierr) 
      call TaoSetObjectiveAndGradientRoutine(tao,FormFunctionGradient,0,ierr);CHKERRA(ierr) 
      call TaoSetObjectiveRoutine(tao,FormFunction,0,ierr);CHKERRA(ierr) 
      call TaoSetGradientRoutine(tao,FormGradient,0,ierr);CHKERRA(ierr) 
      if (hessian) then
        call TaoSetHessianRoutine(tao,Hess,Hess,FormHessian,0,ierr);CHKERRA(ierr)
      endif

      call TaoSetFromOptions(tao,ierr);CHKERRA(ierr) 

      ! Set a minimum tolerance
      call xertst (tol_opt.gt.0.0, 'faulty internal parameter tol_opt')
      call TaoSetTolerances(tao,tol_opt,tol_opt,tol_opt,ierr);CHKERRA(ierr)

      ! Set max number of iterations
      call xertst (maxiter.gt.0, 'faulty internal parameter maxiter')
      call TaoSetMaximumIterations(tao, maxiter, ierr);CHKERRA(ierr)
      call TaoSetMaximumFunctionEvaluations(tao, 2*maxiter, ierr);CHKERRA(ierr)

      ! Solve
      call TaoSolve(tao,ierr);CHKERRA(ierr) 
      call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'PETSC-TAO.OUT',viewer,ierr);CHKERRA(ierr)
      call PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_INFO_DETAIL,ierr);CHKERRA(ierr)
      call TaoView(tao,viewer,ierr);CHKERRA(ierr)
      call PetscViewerDestroy(viewer,ierr);CHKERRA(ierr)

      ! Finalize Memory
      call DestroyProblem(ierr);CHKERRA(ierr)

      call TaoDestroy(tao,ierr);CHKERRA(ierr)
      call PetscFinalize(ierr)

      call b2mn_fin_d
      deallocate(par_opt_phys)
#ifdef TGT
      deallocate(par_opt_physd)
#endif
#ifdef ADJ
      deallocate(par_opt_physb)
#endif
      stop 'b2optim'

      stop

      contains

      subroutine InitializeProblem(npar,ierr)
      implicit none
      PetscReal zero
      PetscErrorCode ierr
      integer ipar, npar

      zero = 0.0

      call VecCreateSeq(PETSC_COMM_SELF,npar,X,ierr);CHKERRQ(ierr)
      call VecDuplicate(X,X_L,ierr);CHKERRQ(ierr)
      call VecDuplicate(X,X_U,ierr);CHKERRQ(ierr)
      call VecSet(X,zero,ierr);CHKERRQ(ierr)
      call VecSet(X_L,zero,ierr);CHKERRQ(ierr)
      call VecSet(X_U,zero,ierr);CHKERRQ(ierr)

      do ipar = 1, npar
        call VecSetValue(X,ipar-1,x0(ipar),INSERT_VALUES,ierr);CHKERRQ(ierr)
        if (xl(ipar).lt.-inf_opt) then
          write(*,*) 'TAO: warning, X_L(',ipar,') set to infty'
          call VecSetValue(X_L,ipar-1,PETSC_NINFINITY,INSERT_VALUES,ierr);CHKERRQ(ierr)
        else
          call VecSetValue(X_L,ipar-1,xl(ipar),INSERT_VALUES,ierr);CHKERRQ(ierr)
        endif
        if (xu(ipar).gt.inf_opt) then
          write(*,*) 'TAO: warning, X_U(',ipar,') set to infty'
          call VecSetValue(X_U,ipar-1,PETSC_INFINITY,INSERT_VALUES,ierr);CHKERRQ(ierr)
        else
          call VecSetValue(X_U,ipar-1,xu(ipar),INSERT_VALUES,ierr);CHKERRQ(ierr)
        endif
      end do
      call VecAssemblyBegin(x,ierr);CHKERRQ(ierr)
      call VecAssemblyEnd(x,ierr);CHKERRQ(ierr)
      write(*,*) 'TAO SET X0'
      call VecView(X,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
      write(*,*) 'TAO SET X_L'
      call VecView(X_L,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
      write(*,*) 'TAO SET X_U'
      call VecView(X_U,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
      ierr = 0

!  Allocate storage space for Hessian;
      call MatCreateSeqBAIJ(PETSC_COMM_SELF,npar,npar,npar,1,PETSC_NULL_INTEGER,Hess,ierr);CHKERRQ(ierr)
      call MatSetOption(Hess,MAT_SYMMETRIC,PETSC_TRUE,ierr);CHKERRQ(ierr)

      end subroutine InitializeProblem


      subroutine DestroyProblem(ierr)
      use taomodule ! IGNORE
      implicit none

      PetscErrorCode ierr

      call VecDestroy(X,ierr);CHKERRQ(ierr)
      call VecDestroy(X_L,ierr);CHKERRQ(ierr)
      call VecDestroy(X_U,ierr);CHKERRQ(ierr)
      call MatDestroy(Hess,ierr);CHKERRQ(ierr)
      ierr = 0
      end subroutine DestroyProblem

      subroutine FormFunctionGradient(tao, XX, F, grad, dummy, ierr)
      implicit none
      real(kind=r8) j(nncf), jd(nncf), gradd(npar_opt)
      integer ipar, isigma
      character*3 str
      PetscErrorCode ierr
      PetscInt dummy
      Vec XX,grad
      Tao tao
      PetscScalar F
      PetscScalar x_v(0:npar_opt-1),g_v(0:npar_opt-1)
      PetscOffset x_i,g_i


      call VecGetArrayRead(XX,x_v,x_i,ierr);CHKERRQ(ierr)
      call VecGetArray(grad,g_v,g_i,ierr);CHKERRQ(ierr)
!     f=(x1-2)^2 + (x2-2)^2 -2*x1-2*x2
!     f=(x_v(x_i)-2.0)*(x_v(x_i)-2.0)+(x_v(x_i+1)-2.0)*(x_v(x_i+1)-2.0)-2.0*(x_v(x_i)+x_v(x_i+1))
!     gx1=2*(x1-2) -2
!     g_v(g_i) = 2.0*(x_v(x_i)-2.0) - 2.0
!     gx2=2*(x2-2) -2
!     g_v(g_i+1) = 2.0*(x_v(x_i+1)-2.0) - 2.0
      do ipar = 1, npar_opt - nsigma_opt
        par_opt_phys(ipar) = x_v(x_i+ipar-1)
        write(str,"(I1)") ipar
        if (ipar.ge.10) write(str,"(I2)") ipar
        write(*,*) 'TAO: eval_F_grad_F with x',trim(str),'= ', par_opt_phys(ipar)
      end do
      isigma = npar_opt - nsigma_opt + 1
      do ipar = 1, nsigma
        if (sigma_opt(ipar)) then
          sigma(ipar) = x_v(x_i+isigma-1)
          write(str,"(I1)") isigma
          if (isigma.ge.10) write(str,"(I2)") isigma
          write(*,*) 'TAO: eval_F_grad_F with x',trim(str),'= ', sigma(ipar)
          isigma = isigma + 1
        endif
      end do
      call b2mn_step_d(j,jd,switchdiff)
      F = j(1)
#ifdef TGT
      do ipar = 1, npar_opt
        g_v(g_i+ipar-1) = jd(1)/par_rescale(ipar) ! rescale par to get order unity
      end do
#endif
#ifdef ADJ
      call set_adj_gradient(npar_opt,gradd,switchdiff)
      do ipar = 1, npar_opt
        g_v(g_i+ipar-1) = gradd(ipar)/par_rescale(ipar) ! rescale par to get order unity
      end do
#endif
      call VecRestoreArrayRead(XX,x_v,x_i,ierr);CHKERRQ(ierr)
      call VecRestoreArray(grad,g_v,g_i,ierr);CHKERRQ(ierr)
      ierr = 0
      end subroutine FormFunctionGradient

      subroutine FormFunction(tao, XX, F, dummy, ierr)
      implicit none
      real(kind=r8) j(nncf)
      integer ipar, isigma
      character*3 str
      PetscErrorCode ierr
      PetscInt dummy
      Vec XX
      Tao tao
      PetscScalar F
      PetscScalar x_v(0:npar_opt-1)
      PetscOffset x_i


      call VecGetArrayRead(XX,x_v,x_i,ierr);CHKERRQ(ierr)

      do ipar = 1, npar_opt - nsigma_opt
        par_opt_phys(ipar) = x_v(x_i+ipar-1)
        write(str,"(I1)") ipar
        if (ipar.ge.10) write(str,"(I2)") ipar
        write(*,*) 'TAO: eval_F with x',trim(str),'= ', par_opt_phys(ipar)
      end do
      isigma = npar_opt - nsigma_opt + 1
      do ipar = 1, nsigma
        if (sigma_opt(ipar)) then
          sigma(ipar) = x_v(x_i+isigma-1)
          write(str,"(I1)") isigma
          if (isigma.ge.10) write(str,"(I2)") isigma
          write(*,*) 'TAO: eval_F with x',trim(str),'= ', sigma(ipar)
          isigma = isigma + 1
        endif
      end do
      call b2mn_step(j)
      F = j(1)

      call VecRestoreArrayRead(XX,x_v,x_i,ierr);CHKERRQ(ierr)
      ierr = 0
      end subroutine FormFunction

      subroutine FormGradient(tao, XX, grad, dummy, ierr)
      implicit none
      real(kind=r8) j(nncf), jd(nncf), gradd(npar_opt)
      integer ipar, isigma
      character*3 str
      PetscErrorCode ierr
      PetscInt dummy
      Vec XX,grad
      Tao tao
      PetscScalar x_v(0:npar_opt-1),g_v(0:npar_opt-1)
      PetscOffset x_i,g_i


      call VecGetArrayRead(XX,x_v,x_i,ierr);CHKERRQ(ierr)
      call VecGetArray(grad,g_v,g_i,ierr);CHKERRQ(ierr)

      do ipar = 1, npar_opt - nsigma_opt
        par_opt_phys(ipar) = x_v(x_i+ipar-1)
        write(str,"(I1)") ipar
        if (ipar.ge.10) write(str,"(I2)") ipar
        write(*,*) 'TAO: eval_grad_F with x',trim(str),'= ', par_opt_phys(ipar)
      end do
      isigma = npar_opt - nsigma_opt + 1
      do ipar = 1, nsigma
        if (sigma_opt(ipar)) then
          sigma(ipar) = x_v(x_i+isigma-1)
          write(str,"(I1)") isigma
          if (isigma.ge.10) write(str,"(I2)") isigma
          write(*,*) 'TAO: eval_grad_F with x',trim(str),'= ', sigma(ipar)
          isigma = isigma + 1
        endif
      end do
      call b2mn_step_d(j,jd,switchdiff)
#ifdef TGT
      do ipar = 1, npar_opt
        g_v(g_i+ipar-1) = jd(1)/par_rescale(ipar) ! rescale par to get order unity
      end do
#endif
#ifdef ADJ
      call set_adj_gradient(npar_opt,gradd,switchdiff)
      do ipar = 1, npar_opt
        g_v(g_i+ipar-1) = gradd(ipar)/par_rescale(ipar) ! rescale par to get order unity
      end do
#endif
      call VecRestoreArrayRead(XX,x_v,x_i,ierr);CHKERRQ(ierr)
      call VecRestoreArray(grad,g_v,g_i,ierr);CHKERRQ(ierr)
      ierr = 0
      end subroutine FormGradient


      subroutine FormHessian(tao,XX,HH,PrecH,dummy,ierr)
      implicit none
      Tao              tao
      Vec              XX
      Mat              HH
      Mat              PrecH
      PetscErrorCode   ierr
      PetscInt         dummy
      PetscBool assembled
      PetscInt         i
      PetscReal        v

 !  Zero existing matrix entries
      call MatAssembled(HH,assembled,ierr);CHKERRQ(ierr)
      if (assembled .eqv. PETSC_TRUE) call MatZeroEntries(HH,ierr);CHKERRQ(ierr)

 !  Compute Hessian entries
      do i=0,npar_opt-1
        v = 1.0
        call MatSetValues(HH,1,i,1,i,v,INSERT_VALUES,ierr);CHKERRQ(ierr)
      enddo

 !  Assemble matrix
      call MatAssemblyBegin(HH,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      call MatAssemblyEnd(HH,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      write(*,*) 'TAO SET H'
      call MatView(Hess,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)

      return
      end subroutine FormHessian

      end program b2optim_petsc
