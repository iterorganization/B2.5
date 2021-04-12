      module taomodule
#include <petsc/finclude/petsctao.h>
      use petsctao ! IGNORE

      Vec X,X_L,X_U
      end module taomodule

      program b2optim_petsc
      use taomodule ! IGNORE
      use b2mod_par_opt_diff
      use b2mod_main_diff &
      , only : b2mn_init_d, b2mn_fin_d
      implicit none

      PetscErrorCode       ierr
      Tao                  tao
      KSP                  ksp
      PC                   pc
      PetscReal            tol
      PetscViewer viewer

      integer :: ncon, nele_jac, ipar

      external FormFunctionGradient

      call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
      if (ierr .ne. 0) then
         print*,'Unable to initialize PETSc'
         stop
      endif

      call read_b2mod_par_opt(ncon, nele_jac)

!      call PetscPrintf(PETSC_COMM_SELF,'Solution should be f(1,1)=-2\n',ierr);CHKERRA(ierr)

      call InitializeProblem(npar_opt,ierr);CHKERRA(ierr)

      ! Allocate and initialize par_opt variables to be used in B2.5
      call b2mn_init_d
      allocate(par_opt_phys(npar_opt))
      if (nsigma.gt.0) then
        allocate(par_opt_physd(nbdirsmax, npar_opt))
        par_opt_physd =0.0_R8
        par_opt_physd(2,2) = 1.0_R8
      endif
      par_opt_phys = 0.0_R8
      flag_optim  = .true. !csc this will tell in b2tqna to use par_opt_phys for parm_dna

      call TaoCreate(PETSC_COMM_SELF,tao,ierr);CHKERRA(ierr) 
      call TaoSetType(tao,TAOBQNLS,ierr);CHKERRA(ierr) 
      call TaoSetInitialVector(tao,X,ierr);CHKERRA(ierr)
      call TaoSetVariableBounds(tao,X_L,X_U,ierr);CHKERRA(ierr) 
      call TaoSetObjectiveAndGradientRoutine(tao,FormFunctionGradient,0,ierr);CHKERRA(ierr) 

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
      if (nsigma.gt.0) deallocate(par_opt_physd)
      stop 'b2optim'

      stop
      end program b2optim_petsc


      subroutine InitializeProblem(npar,ierr)
      use taomodule ! IGNORE
      use b2mod_par_opt_diff
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
      end subroutine InitializeProblem


      subroutine DestroyProblem(ierr)
      use taomodule ! IGNORE
      implicit none

      PetscErrorCode ierr

      call VecDestroy(X,ierr);CHKERRQ(ierr)
      call VecDestroy(X_L,ierr);CHKERRQ(ierr)
      call VecDestroy(X_U,ierr);CHKERRQ(ierr)
      ierr = 0
      end subroutine DestroyProblem

      subroutine FormFunctionGradient(tao, XX, F, grad, dummy, ierr)
      use taomodule ! IGNORE
      use b2mod_types
      use b2mod_main_diff &
      , only : b2mn_step_d
      use b2mod_ad_diff &
      , only : nncf
      use b2mod_par_opt_diff
      implicit none
      real(kind=r8) j(nncf), jd(nncf)
      integer ipar
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
      do ipar = 1, npar_opt
        par_opt_phys(ipar) = x_v(x_i+ipar-1)
        write(str,"(I1)") ipar
        if (ipar.ge.10) write(str,"(I2)") ipar
        write(*,*) 'TAO: eval_F_grad_F with x',trim(str),'= ', par_opt_phys(ipar)
      end do
      call b2mn_step_d(j,jd)
      F = j(1)
      do ipar = 1, npar_opt
        g_v(g_i+ipar-1) = jd(1)
      end do

      call VecRestoreArrayRead(XX,x_v,x_i,ierr);CHKERRQ(ierr)
      call VecRestoreArray(grad,g_v,g_i,ierr);CHKERRQ(ierr)
      ierr = 0
      end subroutine FormFunctionGradient
