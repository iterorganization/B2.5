      module b2mod_par_opt
      use b2mod_types
      implicit none
      integer, save :: npar_opt
      real (kind=R8), allocatable, save :: par_opt0(:), par_opt(:),
     &       par_opt_phys(:)
      logical, save :: flag_optim = .false.
      integer nvmx, nnvar, nncon, nnjac, nsigmx
      parameter (nvmx = 20, nsigmx=10)
      double precision :: x0(nvmx), xl(nvmx), xu(nvmx), gl(nvmx),
     &  gu(nvmx), jcol(nvmx*nvmx), jrow(nvmx*nvmx), jj(nvmx*nvmx)
! sc some variables that allows to change optimization parameters, read from b2.optimization.parameters
      character*256, save :: filename, limited_memory_update_type,
     &  hessian_approximation 
      double precision, save :: cpu_opt=0.0, tol_opt=1.0e-7
      double precision, parameter :: inf_opt = 1.0e30
      integer, save :: maxiter=100, partype(nvmx), nsigma
      integer, save :: prior_type(nvmx) !sc type of prior
      real (kind=r8), save :: prior_par(nvmx,2), prior_range(nvmx,2) !parameters for prior
      logical :: file_ok
      private optimization
      namelist /optimization/ nnvar, nncon, nnjac, x0, xl, xu, gl, gu, 
     & jcol, jrow, jj, maxiter, cpu_opt, tol_opt, hessian_approximation,
     & limited_memory_update_type, partype, prior_type, prior_par,
     & prior_range, nsigma

      contains

      subroutine read_b2mod_par_opt(ncon, nele_jac)
      use b2mod_types
      implicit none
      integer ncon, nele_jac, ii

      filename='b2.optimization.parameters'
      write (*,*) 'OPTIM: max number of readable parameters '//
     & 'and constraints:',nvmx
      write (*,*) 'OPTIM: infinity bound set to:',inf_opt
      x0 = inf_opt*10.0
      xl = -inf_opt*10.0
      xu = inf_opt*10.0
      gl = -inf_opt*10.0
      gu = inf_opt*10.0
      NELE_JAC = 0
      jcol = 0
      jrow = 0
      jj = 0.0
      partype = -1
      nsigma = 0
      prior_type = -1
      prior_par = -1.0_R8
! sc  Set prior range to (-inf,+inf)
      prior_range(:,1) = -inf_opt*10.0
      prior_range(:,2) = inf_opt*10.0
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
      call xertst(nnvar.le.nvmx,'b2mod_par_opt: increase size of nvmx')
      call xertst(nncon.le.nvmx,'b2mod_par_opt: increase size of nvmx')
      call xertst(nnjac.le.nvmx,'b2mod_par_opt: increase size of nvmx')
      call xertst(nnvar.gt.0,'b2mod_par_opt: nnvar<=0')
      call xertst(nncon.gt.0,'b2mod_par_opt: nncon<=0')
      call xertst(nnjac.gt.0,'b2mod_par_opt: nnjac<=0')
! sc  partype refers to the type of design variable:
!     1 - D, particle diffusion coefficient
!     2 - X_e, electron heat diffusion coefficient
!     3 - X_i, ion heat diffusion coefficient
!     4 - sigma, standard deviation on parameters for Bayesian inference. They must be at the last parameters in the vector!
      call xertst(any(partype(1:nnvar).gt.0).and.
     &  any(partype(1:nnvar).le.4),
     &  'b2mod_par_opt: partype<=0 or partype>=4')
! sc  Suppose nsigma switches on Bayesian MAP
      if (nsigma.gt.0) then
! sc check if sigmas are at the end of the parameter vector, if not, issue error
        do ii = nnvar-nsigma+1, nnvar
          call xertst(partype(ii).eq.4, 'b2mod_par_opt: partype 4 '//
     &     'variables MUST be at the end of the vector')
        end do
! sc  prior_type refers to the prior type for the design variables:
!     1 - Uninformative, proper, Gaussian prior.
!     2 - Uniform distribution.
        call xertst(any(prior_type(1:nnvar).gt.0).and.
     &   any(prior_type(1:nnvar).le.2),
     &   'b2mod_par_opt: prior_type<=0 or prior_type>=2')
! sc  prior_par refers to the prior parameters:
!     For prior_type = 1, prior_par(ii,1) = mu, prior_par(ii,2) = sigma
!     For prior_type = 2, not used?
        do ii = 1, nnvar
          if (prior_type(ii) .eq. 1) 
     &      call xertst(prior_par(ii,2) .gt. 0.0_R8,
     &     'b2mod_par_opt: sigma prior <=0')
          call xertst(prior_range(ii,1).lt.prior_range(ii,2),
     &     'b2mod_par_opt: prior_range(.,1)>prior_range(.,2)')
        end do
      end if
      call xertst(any(jcol(1:nnjac).gt.0),'b2mod_par_opt: jcol<=0')
      call xertst(any(jrow(1:nnjac).gt.0),'b2mod_par_opt: jrow<=0')
      call xertst(maxval(jcol(1:nnjac)).le.nncon,
     &  'b2mod_par_opt: jcol>ncon')
      call xertst(maxval(jrow(1:nnjac)).le.nnvar,
     &  'b2mod_par_opt: jrow>nvar')

      npar_opt = nnvar
      ncon = nncon
      nele_jac = nnjac

      return
      end subroutine read_b2mod_par_opt

      end module b2mod_par_opt
