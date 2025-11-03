*-----------------------------------------------------------------------
*.specification

      program b2op
      use b2mod_optim
      use b2mod_exp
      use b2mod_exp_file_names
      use b2mod_par_opt

      implicit none
#ifdef USE_MPI
      include 'mpif.h'
#endif
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2OP is the main program. The code solves the problem of 

*-----------------------------------------------------------------------
*.declarations
      integer(2) :: iresult
      integer :: number_exp = 0
      integer :: i1, i2, iy, i, j
      integer :: lw                                                      !sv 17.01.09
      real (kind=R8) :: func                                             !sv 17.01.09
      real (kind=R8), allocatable :: w(:), e(:)                          !sv 17.01.09
      character, save :: opt_method*6 = 'no_opt'                         !sv 14.04.09 {
      integer :: nfm, info, kdif, nf
      real (kind=R8) :: step, value_fun 

      integer, allocatable :: iw(:)
      real (kind=R8) :: tol 

      integer maxfev, mode, mp5n, nfev, nprint
      real (kind=R8) :: epsfcn, factor, ftol, gtol, xtol, zero

      integer :: IPARAM(6), LDFJAC
      real (kind=R8) :: RPARAM(7)
      real (kind=R8), allocatable :: XSCALE(:), FSCALE(:), FVEC(:),
     &     FJAC(:,:)
      integer :: ipath = 1, ierr
      real (kind=R8), allocatable :: hess(:,:), grad(:), eigval(:),
     &     eigvec(:,:), wrk1(:), wrk2(:)

*   ..common blocks

*   ..local variables
      character*12 hlp_frm

*   ..procedures
      external b2obfn, b2obeq, b2oblm
      external prgini, prgend
*   ..initialize input/output unit numbers

*-----------------------------------------------------------------------
*.computation
*
*   ..program start-up calls
      call prgini ('b2op')
*   .. get experimental data
*   .. for inner lower
      fnm = 'exp_'//trim(fnm_ne_inner_lower)//'.dat'
      inquire(file=fnm,exist=l_ne_inner_lower)
      if(l_ne_inner_lower) then
       call get_dim(fnm, ndat_ne_inner_lower_exp)
       write(*,"(/' The length of file ',a,' is ',i4/)")
     *       fnm, ndat_ne_inner_lower_exp
       allocate(cy_ne_inner_lower_exp(1:ndat_ne_inner_lower_exp), 
     &             ne_inner_lower_exp(1:ndat_ne_inner_lower_exp)) 
       call read_data(fnm, ndat_ne_inner_lower_exp,  
     &                        cy_ne_inner_lower_exp, ne_inner_lower_exp)
       number_exp = number_exp + ndat_ne_inner_lower_exp
*   ..write experimental data for control
       fnm = 'exp_'//trim(fnm_ne_inner_lower)//'.prt'
       call write_data(fnm, ndat_ne_inner_lower_exp, 
     &                       cy_ne_inner_lower_exp,  ne_inner_lower_exp)
      endif

      fnm = 'exp_'//trim(fnm_te_inner_lower)//'.dat'
      inquire(file=fnm,exist=l_te_inner_lower)
      if(l_te_inner_lower) then
       call get_dim(fnm, ndat_te_inner_lower_exp)
       write(*,"(/' The length of file ',a,' is ',i4/)")
     &           fnm, ndat_te_inner_lower_exp
       allocate(cy_te_inner_lower_exp(1:ndat_te_inner_lower_exp),
     &             te_inner_lower_exp(1:ndat_te_inner_lower_exp))
       call read_data(fnm, ndat_te_inner_lower_exp, 
     &                 cy_te_inner_lower_exp, te_inner_lower_exp)
       number_exp = number_exp + ndat_te_inner_lower_exp
*   ..write experimental data for control
       fnm = 'exp_'//trim(fnm_te_inner_lower)//'.prt'
       call write_data(fnm, ndat_te_inner_lower_exp,  
     &                   cy_te_inner_lower_exp, te_inner_lower_exp)
      endif

      fnm = 'exp_'//trim(fnm_ti_inner_lower)//'.dat'
      inquire(file=fnm,exist=l_ti_inner_lower)
      if(l_ti_inner_lower) then
       call get_dim(fnm, ndat_ti_inner_lower_exp)
       write(*,"(/' The length of file ',a,' is ',i4/)")
     &          fnm, ndat_ti_inner_lower_exp
       allocate(cy_ti_inner_lower_exp(1:ndat_ti_inner_lower_exp),
     &             ti_inner_lower_exp(1:ndat_ti_inner_lower_exp))
       call read_data(fnm, ndat_ti_inner_lower_exp, 
     &                 cy_ti_inner_lower_exp, ti_inner_lower_exp)
       number_exp = number_exp + ndat_ti_inner_lower_exp
*   ..write experimental data for control
       fnm = 'exp_'//trim(fnm_ti_inner_lower)//'.prt'
       call write_data(fnm, ndat_ti_inner_lower_exp,  
     &                      cy_ti_inner_lower_exp, ti_inner_lower_exp)
      endif

      fnm = 'exp_'//trim(fnm_js_inner_lower)//'.dat'
      inquire(file=fnm,exist=l_js_inner_lower)
      if(l_js_inner_lower) then
       call get_dim(fnm, ndat_js_inner_lower_exp)
       write(*,"(/' The length of file ',a,' is ',i4/)")
     &          fnm, ndat_js_inner_lower_exp
       allocate(cy_js_inner_lower_exp(1:ndat_js_inner_lower_exp),
     &             js_inner_lower_exp(1:ndat_js_inner_lower_exp))
       call read_data(fnm, ndat_js_inner_lower_exp, 
     &                 cy_js_inner_lower_exp, js_inner_lower_exp)
       number_exp = number_exp + ndat_js_inner_lower_exp
*   ..write experimental data for control
       fnm = 'exp_'//trim(fnm_js_inner_lower)//'.prt'
       call write_data(fnm, ndat_js_inner_lower_exp,  
     &                      cy_js_inner_lower_exp, js_inner_lower_exp)
      endif

*   .. for inner upper
      fnm = 'exp_'//trim(fnm_ne_inner_upper)//'.dat'
      inquire(file=fnm,exist=l_ne_inner_upper)
      if(l_ne_inner_upper) then
       call get_dim(fnm, ndat_ne_inner_upper_exp)
       write(*,"(/' The length of file ',a,' is ',i4/)")
     *       fnm, ndat_ne_inner_upper_exp
       allocate(cy_ne_inner_upper_exp(1:ndat_ne_inner_upper_exp), 
     &             ne_inner_upper_exp(1:ndat_ne_inner_upper_exp)) 
       call read_data(fnm, ndat_ne_inner_upper_exp,  
     &                        cy_ne_inner_upper_exp, ne_inner_upper_exp)
       number_exp = number_exp + ndat_ne_inner_upper_exp
*   ..write experimental data for control
       fnm = 'exp_'//trim(fnm_ne_inner_upper)//'.prt'
       call write_data(fnm, ndat_ne_inner_upper_exp, 
     &                       cy_ne_inner_upper_exp,  ne_inner_upper_exp)
      endif

      fnm = 'exp_'//trim(fnm_te_inner_upper)//'.dat'
      inquire(file=fnm,exist=l_te_inner_upper)
      if(l_te_inner_upper) then
       call get_dim(fnm, ndat_te_inner_upper_exp)
       write(*,"(/' The length of file ',a,' is ',i4/)")
     &           fnm, ndat_te_inner_upper_exp
       allocate(cy_te_inner_upper_exp(1:ndat_te_inner_upper_exp),
     &             te_inner_upper_exp(1:ndat_te_inner_upper_exp))
       call read_data(fnm, ndat_te_inner_upper_exp, 
     &                 cy_te_inner_upper_exp, te_inner_upper_exp)
       number_exp = number_exp + ndat_te_inner_upper_exp
*   ..write experimental data for control
       fnm = 'exp_'//trim(fnm_te_inner_upper)//'.prt'
       call write_data(fnm, ndat_te_inner_upper_exp,  
     &                   cy_te_inner_upper_exp, te_inner_upper_exp)
      endif

      fnm = 'exp_'//trim(fnm_ti_inner_upper)//'.dat'
      inquire(file=fnm,exist=l_ti_inner_upper)
      if(l_ti_inner_upper) then
       call get_dim(fnm, ndat_ti_inner_upper_exp)
       write(*,"(/' The length of file ',a,' is ',i4/)")
     &          fnm, ndat_ti_inner_upper_exp
       allocate(cy_ti_inner_upper_exp(1:ndat_ti_inner_upper_exp),
     &             ti_inner_upper_exp(1:ndat_ti_inner_upper_exp))
       call read_data(fnm, ndat_ti_inner_upper_exp, 
     &                 cy_ti_inner_upper_exp, ti_inner_upper_exp)
       number_exp = number_exp + ndat_ti_inner_upper_exp
*   ..write experimental data for control
       fnm = 'exp_'//trim(fnm_ti_inner_upper)//'.prt'
       call write_data(fnm, ndat_ti_inner_upper_exp,  
     &                      cy_ti_inner_upper_exp, ti_inner_upper_exp)
      endif

      fnm = 'exp_'//trim(fnm_js_inner_upper)//'.dat'
      inquire(file=fnm,exist=l_js_inner_upper)
      if(l_js_inner_upper) then
       call get_dim(fnm, ndat_js_inner_upper_exp)
       write(*,"(/' The length of file ',a,' is ',i4/)")
     &          fnm, ndat_js_inner_upper_exp
       allocate(cy_js_inner_upper_exp(1:ndat_js_inner_upper_exp),
     &             js_inner_upper_exp(1:ndat_js_inner_upper_exp))
       call read_data(fnm, ndat_js_inner_upper_exp, 
     &                 cy_js_inner_upper_exp, js_inner_upper_exp)
       number_exp = number_exp + ndat_js_inner_upper_exp
*   ..write experimental data for control
       fnm = 'exp_'//trim(fnm_js_inner_upper)//'.prt'
       call write_data(fnm, ndat_js_inner_upper_exp,  
     &                      cy_js_inner_upper_exp, js_inner_upper_exp)
      endif

*   .. for outer midplane
      fnm = 'exp_'//trim(fnm_ne_outer_mdpl)//'.dat'
      inquire(file=fnm,exist=l_ne_outer_mdpl)
      if(l_ne_outer_mdpl) then
       call get_dim(fnm, ndat_ne_outer_mdpl_exp)
       write(*,"(/' The length of file ',a,' is ',i4/)")
     *       fnm, ndat_ne_outer_mdpl_exp
       allocate(cy_ne_outer_mdpl_exp(1:ndat_ne_outer_mdpl_exp), 
     &             ne_outer_mdpl_exp(1:ndat_ne_outer_mdpl_exp)) 
       call read_data(fnm, ndat_ne_outer_mdpl_exp, cy_ne_outer_mdpl_exp, 
     &                                 ne_outer_mdpl_exp)
       number_exp = number_exp + ndat_ne_outer_mdpl_exp
*   ..write experimental data for control
       fnm = 'exp_'//trim(fnm_ne_outer_mdpl)//'.prt'
       call write_data(fnm, ndat_ne_outer_mdpl_exp, 
     &                        cy_ne_outer_mdpl_exp,  ne_outer_mdpl_exp)
      endif

      fnm = 'exp_'//trim(fnm_te_outer_mdpl)//'.dat'
      inquire(file=fnm,exist=l_te_outer_mdpl)
      if(l_te_outer_mdpl) then
       call get_dim(fnm, ndat_te_outer_mdpl_exp)
       write(*,"(/' The length of file ',a,' is ',i4/)")
     &           fnm, ndat_te_outer_mdpl_exp
       allocate(cy_te_outer_mdpl_exp(1:ndat_te_outer_mdpl_exp),
     &             te_outer_mdpl_exp(1:ndat_te_outer_mdpl_exp))
       call read_data(fnm, ndat_te_outer_mdpl_exp, cy_te_outer_mdpl_exp,
     &                                 te_outer_mdpl_exp)
       number_exp = number_exp + ndat_te_outer_mdpl_exp
*   ..write experimental data for control
       fnm = 'exp_'//trim(fnm_te_outer_mdpl)//'.prt'
       call write_data(fnm, ndat_te_outer_mdpl_exp,  
     &                   cy_te_outer_mdpl_exp, te_outer_mdpl_exp)
      endif

      fnm = 'exp_'//trim(fnm_ti_outer_mdpl)//'.dat'
      inquire(file=fnm,exist=l_ti_outer_mdpl)
      if(l_ti_outer_mdpl) then
       call get_dim(fnm, ndat_ti_outer_mdpl_exp)
       write(*,"(/' The length of file ',a,' is ',i4/)")
     &          fnm, ndat_ti_outer_mdpl_exp
       allocate(cy_ti_outer_mdpl_exp(1:ndat_ti_outer_mdpl_exp),
     &             ti_outer_mdpl_exp(1:ndat_ti_outer_mdpl_exp))
       call read_data(fnm, ndat_ti_outer_mdpl_exp, cy_ti_outer_mdpl_exp,
     &                                 ti_outer_mdpl_exp)
       number_exp = number_exp + ndat_ti_outer_mdpl_exp
*   ..write experimental data for control
       fnm = 'exp_'//trim(fnm_ti_outer_mdpl)//'.prt'
       call write_data(fnm, ndat_ti_outer_mdpl_exp,  
     &                      cy_ti_outer_mdpl_exp, ti_outer_mdpl_exp)
      endif

*   .. for outer lower
      fnm = 'exp_'//trim(fnm_ne_outer_lower)//'.dat'
      inquire(file=fnm,exist=l_ne_outer_lower)
      if(l_ne_outer_lower) then
       call get_dim(fnm, ndat_ne_outer_lower_exp)
       write(*,"(/' The length of file ',a,' is ',i4/)")
     *       fnm, ndat_ne_outer_lower_exp
       allocate(cy_ne_outer_lower_exp(1:ndat_ne_outer_lower_exp), 
     &             ne_outer_lower_exp(1:ndat_ne_outer_lower_exp)) 
       call read_data(fnm, ndat_ne_outer_lower_exp,  
     &                        cy_ne_outer_lower_exp, ne_outer_lower_exp)
       number_exp = number_exp + ndat_ne_outer_lower_exp
*   ..write experimental data for control
       fnm = 'exp_'//trim(fnm_ne_outer_lower)//'.prt'
       call write_data(fnm, ndat_ne_outer_lower_exp, 
     &                       cy_ne_outer_lower_exp,  ne_outer_lower_exp)
      endif

      fnm = 'exp_'//trim(fnm_te_outer_lower)//'.dat'
      inquire(file=fnm,exist=l_te_outer_lower)
      if(l_te_outer_lower) then
       call get_dim(fnm, ndat_te_outer_lower_exp)
       write(*,"(/' The length of file ',a,' is ',i4/)")
     &           fnm, ndat_te_outer_lower_exp
       allocate(cy_te_outer_lower_exp(1:ndat_te_outer_lower_exp),
     &             te_outer_lower_exp(1:ndat_te_outer_lower_exp))
       call read_data(fnm, ndat_te_outer_lower_exp, 
     &                 cy_te_outer_lower_exp, te_outer_lower_exp)
       number_exp = number_exp + ndat_te_outer_lower_exp
*   ..write experimental data for control
       fnm = 'exp_'//trim(fnm_te_outer_lower)//'.prt'
       call write_data(fnm, ndat_te_outer_lower_exp,  
     &                   cy_te_outer_lower_exp, te_outer_lower_exp)
      endif

      fnm = 'exp_'//trim(fnm_ti_outer_lower)//'.dat'
      inquire(file=fnm,exist=l_ti_outer_lower)
      if(l_ti_outer_lower) then
       call get_dim(fnm, ndat_ti_outer_lower_exp)
       write(*,"(/' The length of file ',a,' is ',i4/)")
     &          fnm, ndat_ti_outer_lower_exp
       allocate(cy_ti_outer_lower_exp(1:ndat_ti_outer_lower_exp),
     &             ti_outer_lower_exp(1:ndat_ti_outer_lower_exp))
       call read_data(fnm, ndat_ti_outer_lower_exp, 
     &                 cy_ti_outer_lower_exp, ti_outer_lower_exp)
       number_exp = number_exp + ndat_ti_outer_lower_exp
*   ..write experimental data for control
       fnm = 'exp_'//trim(fnm_ti_outer_lower)//'.prt'
       call write_data(fnm, ndat_ti_outer_lower_exp,  
     &                      cy_ti_outer_lower_exp, ti_outer_lower_exp)
      endif

      fnm = 'exp_'//trim(fnm_js_outer_lower)//'.dat'
      inquire(file=fnm,exist=l_js_outer_lower)
      if(l_js_outer_lower) then
       call get_dim(fnm, ndat_js_outer_lower_exp)
       write(*,"(/' The length of file ',a,' is ',i4/)")
     &          fnm, ndat_js_outer_lower_exp
       allocate(cy_js_outer_lower_exp(1:ndat_js_outer_lower_exp),
     &             js_outer_lower_exp(1:ndat_js_outer_lower_exp))
       call read_data(fnm, ndat_js_outer_lower_exp, 
     &                 cy_js_outer_lower_exp, js_outer_lower_exp)
       number_exp = number_exp + ndat_js_outer_lower_exp
*   ..write experimental data for control
       fnm = 'exp_'//trim(fnm_js_outer_lower)//'.prt'
       call write_data(fnm, ndat_js_outer_lower_exp,  
     &                      cy_js_outer_lower_exp, js_outer_lower_exp)
      endif

*   .. for outer upper
      fnm = 'exp_'//trim(fnm_ne_outer_upper)//'.dat'
      inquire(file=fnm,exist=l_ne_outer_upper)
      if(l_ne_outer_upper) then
       call get_dim(fnm, ndat_ne_outer_upper_exp)
       write(*,"(/' The length of file ',a,' is ',i4/)")
     *       fnm, ndat_ne_outer_upper_exp
       allocate(cy_ne_outer_upper_exp(1:ndat_ne_outer_upper_exp), 
     &             ne_outer_upper_exp(1:ndat_ne_outer_upper_exp)) 
       call read_data(fnm, ndat_ne_outer_upper_exp,  
     &                        cy_ne_outer_upper_exp, ne_outer_upper_exp)
       number_exp = number_exp + ndat_ne_outer_upper_exp
*   ..write experimental data for control
       fnm = 'exp_'//trim(fnm_ne_outer_upper)//'.prt'
       call write_data(fnm, ndat_ne_outer_upper_exp, 
     &                       cy_ne_outer_upper_exp,  ne_outer_upper_exp)
      endif

      fnm = 'exp_'//trim(fnm_te_outer_upper)//'.dat'
      inquire(file=fnm,exist=l_te_outer_upper)
      if(l_te_outer_upper) then
       call get_dim(fnm, ndat_te_outer_upper_exp)
       write(*,"(/' The length of file ',a,' is ',i4/)")
     &           fnm, ndat_te_outer_upper_exp
       allocate(cy_te_outer_upper_exp(1:ndat_te_outer_upper_exp),
     &             te_outer_upper_exp(1:ndat_te_outer_upper_exp))
       call read_data(fnm, ndat_te_outer_upper_exp, 
     &                 cy_te_outer_upper_exp, te_outer_upper_exp)
       number_exp = number_exp + ndat_te_outer_upper_exp
*   ..write experimental data for control
       fnm = 'exp_'//trim(fnm_te_outer_upper)//'.prt'
       call write_data(fnm, ndat_te_outer_upper_exp,  
     &                   cy_te_outer_upper_exp, te_outer_upper_exp)
      endif

      fnm = 'exp_'//trim(fnm_ti_outer_upper)//'.dat'
      inquire(file=fnm,exist=l_ti_outer_upper)
      if(l_ti_outer_upper) then
       call get_dim(fnm, ndat_ti_outer_upper_exp)
       write(*,"(/' The length of file ',a,' is ',i4/)")
     &          fnm, ndat_ti_outer_upper_exp
       allocate(cy_ti_outer_upper_exp(1:ndat_ti_outer_upper_exp),
     &             ti_outer_upper_exp(1:ndat_ti_outer_upper_exp))
       call read_data(fnm, ndat_ti_outer_upper_exp, 
     &                 cy_ti_outer_upper_exp, ti_outer_upper_exp)
       number_exp = number_exp + ndat_ti_outer_upper_exp
*   ..write experimental data for control
       fnm = 'exp_'//trim(fnm_ti_outer_upper)//'.prt'
       call write_data(fnm, ndat_ti_outer_upper_exp,  
     &                      cy_ti_outer_upper_exp, ti_outer_upper_exp)
      endif

      fnm = 'exp_'//trim(fnm_js_outer_upper)//'.dat'
      inquire(file=fnm,exist=l_js_outer_upper)
      if(l_js_outer_upper) then
       call get_dim(fnm, ndat_js_outer_upper_exp)
       write(*,"(/' The length of file ',a,' is ',i4/)")
     &          fnm, ndat_js_outer_upper_exp
       allocate(cy_js_outer_upper_exp(1:ndat_js_outer_upper_exp),
     &             js_outer_upper_exp(1:ndat_js_outer_upper_exp))
       call read_data(fnm, ndat_js_outer_upper_exp, 
     &                 cy_js_outer_upper_exp, js_outer_upper_exp)
       number_exp = number_exp + ndat_js_outer_upper_exp
*   ..write experimental data for control
       fnm = 'exp_'//trim(fnm_js_outer_upper)//'.prt'
       call write_data(fnm, ndat_js_outer_upper_exp,  
     &                      cy_js_outer_upper_exp, js_outer_upper_exp)
      endif

      fnm = 'exp_'//trim(fnm_power)//'.dat'
      inquire(file=fnm,exist=l_power)
      if(l_power) then
       open(70,file=fnm)
       read(70,*) power_exp
       close(70)
       number_exp = number_exp + 1
*   ..write experimental data for control
       open(70,file='exp_'//trim(fnm_power)//'.prt')
       write(70,"(1x,2(es14.6))") power_exp
       close(70)
      endif

      open(71,file='b2optm.prt',access='append')
      nfunc = 0

!      call ipgetc ('b2op', opt_method)
      open(70,file='b2op.dat')
      read(70,'(i3)') npar_opt
      allocate(par_opt0(1:npar_opt), par_opt(1:npar_opt),
     &         par_opt_phys(1:npar_opt))
      read(70,'(a)') opt_method
      read(70,'(i3)') istyle
      read(70,*) tol
      read(70, *) (par_opt0(i1),i1=1,npar_opt)
      close(70)

      write(71,"(' number_exp point =', i6)") number_exp
      write(71,"(' npar_opt         =', i3)") npar_opt
      write(71,'(a)') opt_method
      write(71,"(' istyle for restriction', i3)") istyle
      write(71,"('Tolerance =',es10.2)") tol
      write(71,"(' par_opt0')")
      write(71,'(es14.5)') par_opt0

      select case (opt_method)

      case('UNCMND')
      write(71,"(' UNCMND')")
      lw=npar_opt*(npar_opt+10)
      allocate(w(lw))
      CALL UNCMND (npar_opt,par_opt0,b2obfn,par_opt,func,INFO,W,LW)
      deallocate(w)
      if(info.eq.0) write(71,"(a)") 'INFO = 0: Optimal solution found'            
      if(info.eq.1) write(71,"(2a)") 'INFO = 1: Terminated with ',
     &           'gradient small, X is probably optimal'           
      if(info.eq.2) write(71,"(2a)") 'INFO = 2: Terminated with ',
     &            'stepsize small, X is probably optimal'          
      if(info.eq.3)write(71,"(2a)") 'INFO = 3: Lower point cannot be',
     &            ' found, X is probably optimal'            
      if(info.eq.4)write(71,"(2a)") 'INFO = 4: Iteration limit (150) ',
     &            'exceeded'            
      if(info.eq.5) write(71,"(2a)") 'INFO = 5: Too many large steps,',
     &             'function may be unbounded'            
      if(info.eq.-1)write(71,"(a)") 'INFO =-1: Insufficient workspace'

      call b2obfn(npar_opt, par_opt, func)

      case('dgz1m')
       write(71,"(' dgz1m')")
       allocate(e(npar_opt))
       info = 1
       nfm  = 200
       step = 0.5_R8
!       tol=1.0e-2_R8
       par_opt = par_opt0
       call dgz1m(npar_opt, par_opt, tol, value_fun, step, e, nfm, 
     &            info, b2obfn)
       write(71,"(' INFO =',i3,' func =',es10.2)") info, value_fun
       call b2obfn(npar_opt, par_opt, func)
       deallocate(e)
#ifdef WINDOWS
      case('DUNLSF')
       write(71,"(' DUNLSF')")
       allocate ( XSCALE(npar_opt), FSCALE(number_exp), FVEC(number_exp),
     *             FJAC(number_exp,npar_opt) )
       LDFJAC = number_exp

       FSCALE = 1.0d0

       CALL DU4LSF (IPARAM, RPARAM)
       IPARAM(2) = 2
       RPARAM(4)=tol    !  1.0e-2_R8

       WRITE (71,*) 'IMSL DU4LSF'


       CALL DUNLSF (b2obeq, Number_exp, npar_opt, par_opt0, XSCALE,
     &     FSCALE, IPARAM, RPARAM, par_opt, FVEC, FJAC, LDFJAC)

       call b2obeq(number_exp, npar_opt, par_opt, fvec)
       deallocate ( XSCALE, FSCALE, FVEC, FJAC )
#endif

      case('lmdif1')
       write(71,"(' lmdif1')")
       lw=number_exp*npar_opt+5*npar_opt+number_exp
       allocate ( iw(npar_opt), w(lw), FVEC(number_exp))

!       tol=1.0e-2_R8
       par_opt = par_opt0
!       call lmdif1(b2oblm, Number_exp, npar_opt, par_opt, fvec, tol,
!     &      info, iw, w, lw)
       factor = 1.0d2
       zero = tol ! 1.0e-2_R8 ! 0.0_R8
       info = 0
       maxfev = 50*(npar_opt + 1) ! 200*(npar_opt + 1)
       ftol = tol
       xtol = tol
       gtol = zero
       epsfcn = 1.0e-2_R8 ! square of numerical difference step
       mode = 1
       nprint = 0
       mp5n = number_exp + 5*npar_opt
       call lmdif(b2oblm, number_exp, npar_opt, par_opt, fvec, ftol, 
     *   xtol, gtol, maxfev, epsfcn,w(1), mode, factor, nprint, info, 
     *   nfev,w(mp5n+1), number_exp, iw, w(npar_opt+1), w(2*npar_opt+1),
     *   w(3*npar_opt+1), w(4*npar_opt+1), w(5*npar_opt+1))
       if (info .eq. 8) info = 4
       if(info.eq.0) write(71,"(2a)")'INFO = 0: ',
     &              'improper input parameters'
       if(info.eq.1) write(71,"(3a)") 'INFO = 1: algorithm estimates ',
     &              'that the relative error between x and',
     &              ' the solution is at most tol'           
       if(info.eq.2) write(71,"(3a)") 'INFO = 2: ',
     &            'algorithm estimates that the relative error',
     &            ' between x and the solution is at most tol'          
       if(info.eq.3)write(71,"(2a)") 'INFO = 3: ',
     &            ' conditions for info = 1 and info = 2 both hold'            
       if(info.eq.4)write(71,"(2a)") 'INFO = 4: fvec is orthogonal to',
     &            ' the columns of the jacobian to machine precision'            
       if(info.eq.5) write(71,"(2a)") 'INFO = 5: number of calls ',
     &             'to fcn has reached or, exceeded 200*(n+1)'            
       if(info.eq.6)write(71,"(2a)") 'INFO =6: tol is too small.',
     & ' no further reduction in the sum of squares is possible'
       if(info.eq.7)write(71,"(2a)") 'INFO =7: tol is too small.',
     & ' no further reduction the approximate solution x is possible'

       call b2oblm(number_exp, npar_opt, par_opt, fvec, i1)

       deallocate (iw, w, FVEC )

      case('dspac8')
       write(71,"(' dspac8')")
       lw=4*npar_opt+2*npar_opt*npar_opt+number_exp*npar_opt+number_exp
       allocate ( iw(npar_opt), w(lw), FVEC(number_exp))
       info = 1
       kdif = 1
       nfm  = 200
       step = 0.5_R8
!       tol=1.0e-2_R8
       par_opt = par_opt0
       call dspac8(npar_opt, number_exp, par_opt, step, tol, 
     *      value_fun, lw, w, iw, nf, nfm, fvec, info, kdif, b2obeq)
       write(71,"(' INFO =',i3,' nf ',i5,' func =',es10.2)") info, nf,
     *         value_fun
       call b2obeq(number_exp, npar_opt, par_opt, fvec)

       deallocate (iw, w, FVEC )

      case('study')
       write(71,"(' Analyze')")
       allocate(hess(npar_opt,npar_opt), grad(npar_opt))
       allocate(eigval(npar_opt), eigvec(npar_opt,npar_opt))
       allocate(wrk1(npar_opt), wrk2(npar_opt))

       call hess_grad(b2obfn, npar_opt, par_opt0, hess, npar_opt, grad)

       write(71,"('    Hesse')")
       write(hlp_frm,'(a,i3,a)') '(',npar_opt,'es12.4)'
       do i =1,npar_opt
        write(71,hlp_frm) (hess(i,j),j=1,npar_opt)
       enddo
       write(71,"('    gradient')")
       do i =1,npar_opt
        write(71,"(5x,es12.4)") grad(i)
       enddo

       call RS(npar_opt, npar_opt, hess, eigval, ipath, eigvec,
     &         wrk1, wrk2, IERR)

       write(71,"(5x,' ierr = ', i5)")ierr
       write(71,"('    eigenvalues ')")
       do i =1,npar_opt
        write(71,"(5x,es16.8)") eigval(i)
       enddo

       write(71,"('    eigenvectors ')")
       write(hlp_frm,'(a,i3,a)') '(',npar_opt,'es16.8)'
       do i =1,npar_opt
        write(71,hlp_frm) (eigvec(i,j),j=1,npar_opt)
       enddo
       open(70,file='Hesse',form='unformatted')
       write(70) hess
       close(70)
       open(70,file='Hesse_evec',form='unformatted')
       write(70) eigvec
       close(70)
       open(70,file='Hesse_eval',form='unformatted')
       write(70) eigval
       close(70)
       if(minval(abs(eigval)).ne.0.0_R8) then
        write(71,"(a,es12.4)")  ' cond(Hess) =',
     &    maxval(abs(eigval))/minval(abs(eigval))

        do i = 1, npar_opt      
         wrk1(i) = 0.0
         do j = 1, npar_opt
          wrk1(i) = wrk1(i) - eigvec(j,i)*grad(j)
         enddo
         wrk1(i)=wrk1(i)/eigval(i)
        enddo

        do i = 1, npar_opt
         wrk2(i) = 0.0
         do j = 1, npar_opt
          wrk2(i) = wrk2(i) + eigvec(i,j)*wrk1(j)
         enddo
         enddo
        write(71,*) 'Delta_par_opt = '
        write(71,"(es15.7)") (wrk2(j), j = 1 , npar_opt)
        par_opt = par_opt0 + wrk2
        write(71,*) 'par_opt = '
        write(71,"(es15.7)") (par_opt(j), j = 1 , npar_opt)
        call b2obfn(npar_opt, par_opt, func)
       endif

       deallocate(hess, grad, eigval, eigvec, wrk1, wrk2)
      case('no_opt')
       write(71,"(' No optimization. ')")
       par_opt = par_opt0

       call b2obfn(npar_opt, par_opt, func)
!       allocate ( FVEC(number_exp))
!!        call b2obeq(number_exp, npar_opt, par_opt, fvec)
!        call b2oblm(number_exp, npar_opt, par_opt, fvec, i1)
!       deallocate ( FVEC)

      case default
       write(71,"(' Do nothing ')")
      end select 
       
      close(71)
      
      deallocate(par_opt0, par_opt, par_opt_phys)
      if(l_ne_outer_mdpl) 
     &   deallocate(cy_ne_outer_mdpl_exp, ne_outer_mdpl_exp) 
      if(l_te_outer_mdpl) 
     &   deallocate(cy_te_outer_mdpl_exp, te_outer_mdpl_exp) 
      if(l_ti_outer_mdpl) 
     &   deallocate(cy_ti_outer_mdpl_exp, ti_outer_mdpl_exp) 

      call prgend ()
      stop 'b2op'

*-----------------------------------------------------------------------
*.end b2op

      end program b2op
