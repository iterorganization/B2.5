      module b2mod_par_opt
      use b2mod_types
      implicit none
      integer :: npar_opt
      real (kind=R8), allocatable, save :: par_opt0(:), par_opt(:),
     &       par_opt_phys(:)
      logical, save :: flag_optim = .false.
      integer nvmx, nnvar, nncon, nnjac
      parameter (nvmx = 10)
      double precision :: x0(nvmx), xl(nvmx), xu(nvmx), gl(nvmx),
     &  gu(nvmx), jcol(nvmx*nvmx), jrow(nvmx*nvmx), jj(nvmx*nvmx)
! sc some variables that allows to change optimization parameters, read from b2.optimization.parameters
      character*256, save :: filename, limited_memory_update_type,
     &  hessian_approximation 
      double precision, save :: cpu_opt=0.0, tol_opt=1.0e-7
      double precision, parameter :: inf_opt = 1.0e30
      integer, save :: maxiter=100
      logical :: file_ok
      private optimization
      namelist /optimization/ nnvar, nncon, nnjac, x0, xl, xu, gl, gu, 
     & jcol, jrow, jj, maxiter, cpu_opt, tol_opt, hessian_approximation,
     & limited_memory_update_type

      contains

      subroutine read_b2mod_par_opt(nvar, ncon, nele_jac)
      use b2mod_types
      implicit none
      integer nvar, ncon, nele_jac

      filename='b2.optimization.parameters'
      write (*,*) 'OPTIM: max number of readablee vars and cons:',nvmx
      write (*,*) 'OPTIM: infinity bound set to:',inf_opt
      x0 = inf_opt*10
      xl = -inf_opt*10
      xu = inf_opt*10
      gl = -inf_opt*10
      gu = inf_opt*10
      NELE_JAC = 0
      jcol = 0
      jrow = 0
      jj = 0.0
      call find_file(filename,file_ok)
      if(file_ok) then
        open(99,file=filename)
        write(*,*) 'Reading namelist optimization from ',trim(filename)
#ifdef USE_MPI
        read(99,optimization,err=999)
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
      call xertst(any(jcol(1:nnjac).gt.0),'b2mod_par_opt: jcol<=0')
      call xertst(any(jrow(1:nnjac).gt.0),'b2mod_par_opt: jrow<=0')
      call xertst(maxval(jcol(1:nnjac)).le.nncon,
     &  'b2mod_par_opt: jcol>ncon')
      call xertst(maxval(jrow(1:nnjac)).le.nnvar,
     &  'b2mod_par_opt: jrow>nvar')

      nvar = nnvar
      ncon = nncon
      nele_jac = nnjac

      return
      end subroutine read_b2mod_par_opt

      end module b2mod_par_opt
