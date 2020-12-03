      module taomodule
#include <petsc/finclude/petsctao.h>
      use petsctao ! IGNORE

      Vec X,X_L,X_U
      end module taomodule

      program b2optim_petsc
      use taomodule ! IGNORE
      use b2mod_par_opt
      use b2mod_main_diff &
      , only : b2mn_init_d, b2mn_fin_d
      implicit none

      PetscErrorCode       ierr
      Tao                  tao
      KSP                  ksp
      PC                   pc
      PetscReal            tol
      PetscViewer viewer

      integer :: nvar, ncon, nele_jac, ivar

      external FormFunctionGradient

      call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
      if (ierr .ne. 0) then
         print*,'Unable to initialize PETSc'
         stop
      endif

      call read_b2mod_par_opt(nvar, ncon, nele_jac)

!      call PetscPrintf(PETSC_COMM_SELF,'Solution should be f(1,1)=-2\n',ierr);CHKERRA(ierr)

      call InitializeProblem(nvar,ierr);CHKERRA(ierr)

      ! Allocate and initialize par_opt variables to be used in B2.5
      call b2mn_init_d
      npar_opt = nvar
      allocate(par_opt_phys(1:npar_opt))
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
      stop 'b2optim'

      stop
      end program b2optim_petsc


      subroutine InitializeProblem(nvar,ierr)
      use taomodule ! IGNORE
      use b2mod_par_opt
      implicit none
      PetscReal zero
      PetscErrorCode ierr
      integer ivar, nvar

      zero = 0.0

      call VecCreateSeq(PETSC_COMM_SELF,nvar,X,ierr);CHKERRQ(ierr)
      call VecDuplicate(X,X_L,ierr);CHKERRQ(ierr)
      call VecDuplicate(X,X_U,ierr);CHKERRQ(ierr)
      call VecSet(X,zero,ierr);CHKERRQ(ierr)
      call VecSet(X_L,zero,ierr);CHKERRQ(ierr)
      call VecSet(X_U,zero,ierr);CHKERRQ(ierr)

      do ivar = 1, nvar
        call VecSetValue(X,ivar-1,x0(ivar),INSERT_VALUES,ierr);CHKERRQ(ierr)
        if (xl(ivar).lt.-inf_opt) then
          write(*,*) 'TAO: warning, X_L(',ivar,') set to infty'
          call VecSetValue(X_L,ivar-1,PETSC_NINFINITY,INSERT_VALUES,ierr);CHKERRQ(ierr)
        else
          call VecSetValue(X_L,ivar-1,xl(ivar),INSERT_VALUES,ierr);CHKERRQ(ierr)
        endif
        if (xu(ivar).gt.inf_opt) then
          write(*,*) 'TAO: warning, X_U(',ivar,') set to infty'
          call VecSetValue(X_U,ivar-1,PETSC_INFINITY,INSERT_VALUES,ierr);CHKERRQ(ierr)
        else
          call VecSetValue(X_U,ivar-1,xu(ivar),INSERT_VALUES,ierr);CHKERRQ(ierr)
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
      use b2mod_user_namelist_diff &
      , only : nncf
      use b2mod_par_opt
      implicit none
      real(kind=r8) j(nncf), jd(nncf)
      PetscErrorCode ierr
      PetscInt dummy
      Vec XX,grad
      Tao tao
      PetscScalar F
      PetscScalar x_v(0:1),g_v(0:1)
      PetscOffset x_i,g_i


      call VecGetArrayRead(XX,x_v,x_i,ierr);CHKERRQ(ierr)
      call VecGetArray(grad,g_v,g_i,ierr);CHKERRQ(ierr)
!     f=(x1-2)^2 + (x2-2)^2 -2*x1-2*x2
!     f=(x_v(x_i)-2.0)*(x_v(x_i)-2.0)+(x_v(x_i+1)-2.0)*(x_v(x_i+1)-2.0)-2.0*(x_v(x_i)+x_v(x_i+1))
!     gx1=2*(x1-2) -2
!     g_v(g_i) = 2.0*(x_v(x_i)-2.0) - 2.0
!     gx2=2*(x2-2) -2
!     g_v(g_i+1) = 2.0*(x_v(x_i+1)-2.0) - 2.0
      par_opt_phys(1) = x_v(x_i) !real(XX(1),kind=r8)
      write(*,*) 'TAO: eval_F_grad_F with x= ', par_opt_phys(1)
      call b2mn_step_d(j,jd)
      F = j(1)
      g_v(g_i) = jd(1) !DBLE(jd(1))

      call VecRestoreArrayRead(XX,x_v,x_i,ierr);CHKERRQ(ierr)
      call VecRestoreArray(grad,g_v,g_i,ierr);CHKERRQ(ierr)
      ierr = 0
      end subroutine FormFunctionGradient
