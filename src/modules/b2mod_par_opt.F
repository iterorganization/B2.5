      module b2mod_par_opt
      use b2mod_types
      use b2mod_user_namelist
     & , only : inf_opt, nsigma, sigma,
     & prior_type, prior_par, prior_range,
     & nsigmx, cftype
      implicit none
      integer, save :: npar_opt = 0
      real (kind=R8), allocatable, save :: par_opt0(:), par_opt(:),
     &       par_opt_phys(:)
      logical, save :: flag_optim = .false.
      integer nvmx, nnvar, nncon, nnjac
      parameter (nvmx = 20)
      double precision :: x0(nvmx), xl(nvmx), xu(nvmx), gl(nvmx),
     &  gu(nvmx), jcol(nvmx*nvmx), jrow(nvmx*nvmx), jj(nvmx*nvmx)
! sc some variables that allows to change optimization parameters, read from b2.optimization.parameters
      character*256, save :: filename, limited_memory_update_type,
     &  hessian_approximation 
      double precision, save :: cpu_opt=0.0, tol_opt=1.0e-7
      integer, save :: maxiter=100, partype(nvmx), nsigma_opt = 0 
      logical, save :: sigma_opt(nsigmx)
      logical :: file_ok
      private optimization
      namelist /optimization/ nnvar, nncon, nnjac, x0, xl, xu, gl, gu, 
     & jcol, jrow, jj, maxiter, cpu_opt, tol_opt, hessian_approximation,
     & limited_memory_update_type, partype, sigma_opt

      contains

      subroutine read_b2mod_par_opt(ncon, nele_jac)
      use b2mod_types
      implicit none
      integer ncon, nele_jac, ii, isigma

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
      nnjac = 0
      jcol = 0
      jrow = 0
      jj = 0.0
      partype = -1
      sigma_opt = .true.
! FIXME set if (read_file)
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
      call xertst(nnvar.gt.0.and.nnvar.le.nvmx,
     & 'b2mod_par_opt: nnvar<=0 or nnvar>nvmx')
      call xertst(nncon.le.nvmx,'b2mod_par_opt: increase size of nvmx')
      call xertst(nnjac.le.nvmx,'b2mod_par_opt: increase size of nvmx')
      call xertst(nnvar.ge.0,'b2mod_par_opt: nnvar<0')
      call xertst(nncon.ge.0,'b2mod_par_opt: nncon<0')
      call xertst(nnjac.ge.0,'b2mod_par_opt: nnjac<0')
!     partype refers to the type of design variable:
!     1 - D, particle diffusion coefficient
!     2 - X_e, electron heat diffusion coefficient
!     3 - X_i, ion heat diffusion coefficient
!     4 - sigma, standard deviation on parameters for Bayesian inference. They must be the last parameters in the vector!
      call xertst(any(partype(1:nnvar).gt.0).and.
     &  any(partype(1:nnvar).le.4),
     &  'b2mod_par_opt: partype<=0 or partype>=4')
      if (nsigma.gt.0 .and. any(sigma_opt(1:nsigma))) then
!      calculate first how many sigmas are being optimized
       nsigma_opt = 0
       do ii = 1, nsigma
        if (sigma_opt(ii)) nsigma_opt = nsigma_opt + 1
       end do
!      check if sigmas are at the end of the parameter vector, if not, issue error
       isigma = nnvar - nsigma_opt + 1 ! position of first sigma in partype
       do ii = 1, nsigma_opt
        call xertst(partype(ii+isigma-1).eq.4,'b2mod_par_opt: '//
     &   'partype 4 variables MUST be at the end of the vector')
       end do
       do ii = 1, nsigma
!       only assign initial value from x0 if that sigma is being optimized
        if (sigma_opt(ii)) then
         sigma(ii) = dble(x0(isigma))
         isigma = isigma + 1 !next sigma in inital vector x0
        endif
       end do
      endif
!     Suppose MAP estimation is turned on only when cftype(1)=6
      if (cftype(1).eq.6) then
!      prior_type refers to the prior type for the design variables:
!      1 - Uninformative, proper, Gaussian prior.
!      2 - Uniform distribution.
       call xertst(any(prior_type(1:nnvar).ge.0).and.
     &  any(prior_type(1:nnvar).le.1),
     &  'b2mod_par_opt: prior_type<=0 or prior_type>=1')
!      prior_par refers to the prior parameters:
!      For prior_type = 1, prior_par(ii,1) = mu, prior_par(ii,2) = sigma
!      For prior_type = 2, not used?
       do ii = 1, nnvar
        if (prior_type(ii) .eq. 1) 
     &   call xertst(prior_par(ii,2) .gt. 0.0_R8,
     &    'b2mod_par_opt: sigma prior <=0')
        call xertst(prior_range(ii,1).lt.prior_range(ii,2),
     &   'b2mod_par_opt: prior_range(.,1)>prior_range(.,2)')
       end do
      end if
      if (nnjac.gt.0) then
       call xertst(any(jcol(1:nnjac).gt.0),'b2mod_par_opt: jcol<=0')
       call xertst(any(jrow(1:nnjac).gt.0),'b2mod_par_opt: jrow<=0')
       call xertst(maxval(jcol(1:nnjac)).le.nncon,
     & 'b2mod_par_opt: jcol>ncon')
       call xertst(maxval(jrow(1:nnjac)).le.nnvar,
     & 'b2mod_par_opt: jrow>nvar')
      endif

      npar_opt = nnvar
      ncon = nncon
      nele_jac = nnjac

      return
      end subroutine read_b2mod_par_opt

      end module b2mod_par_opt
