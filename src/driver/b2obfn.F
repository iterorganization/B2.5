*-----------------------------------------------------------------------
*.specification

      subroutine b2obfn(npar_opt, par_opt, func)
      use b2mod_optim
      use b2mod_exp
      use b2mod_exp_file_names
      use b2mod_par_opt
     , , only : par_opt_phys
#ifdef WINDOWS
      USE DFLIB ! IGNORE
#endif
      implicit none
*   ..input arguments (unchanged on exit)
      integer npar_opt
      real (kind=R8) :: par_opt(1:npar_opt)
*   ..output arguments (unspecified on entry)
      real (kind=R8) :: func

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2OBFN calculates value of objective functional.
*
*
*  2. specification
*
*
*
*
*  3. description (see also routine b2cdca)
*
*
*
*  4. references
*
*
*
*  5. parameters (see also routine b2cdcv)
*
*     None.
*
*
*  6. error indicators
*
*     In case an error condition is detected, a call is made to the
*     routine xerrab. This causes an error message to be printed,
*     after which the program halts.
*
*
*  7. accuracy
*
*     As a measure of convergence we employ the scaled L1-norms of the
*     residuals of the individual equations (continuity for each
*     species, momentum balance for each species, electron heat
*     balance, atom heat balance, and charge balance). The scaling
*     provides numbers that have the dimension of (1/time). A reliable
*     estimate of the accuracy of a calculation can be obtained only
*     by continuing the calculation to better convergence.
*
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks

*   ..local variables
#ifdef WINDOWS
      integer(2) :: iresult
      logical(4) :: lresult
#endif
      logical(4) :: b2fstate_exist
      integer npoint
      real (kind=R8) :: cpu_s, cpu_f, dif, fcain, fab, shift(1:2)
      character*12 hlp_frm
*   ..procedures
      external subini, subend
*-----------------------------------------------------------------------
*.computation
crfs !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
crfs put here the other main programs b2ag (prepares geometry and
crfs  magnetic field), b2ah (prepares tables of default physics
crfs  parameters), b2ai (prepares an initial plasma state),
crfs !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

*   ..program start-up calls
      call subini ('b2obfn')

*   ..calculate physical values of the optimization parameters
!      par_opt_phys(1) = 0.3_R8 + par_opt(1)**2     ! hce
!      par_opt_phys(2) = 0.3_R8 + par_opt(2)**2     ! hci
!      par_opt_phys(3) = 0.3_R8 + par_opt(3)**2     ! Di
      par_opt_phys(1) = fcain(0.3_R8,par_opt(1),istyle)   ! hce
      par_opt_phys(2) = fcain(0.3_R8,par_opt(2),istyle)   ! hci
      par_opt_phys(3) = fcain(0.3_R8,par_opt(3),istyle)   ! Di
!    parameters for boundary conditions
      par_opt_phys(4) = fab(-0.10_R8, -0.001_R8, par_opt(4), istyle)  ! conpar(0,ib,1) [ -0.10,  -0.001]
      par_opt_phys(5) = fab(-0.01_R8, -0.001_R8, par_opt(5), istyle)  ! conpar(1,ib,1) [ -0.01,  -0.001]
      par_opt_phys(6) = fab(-1.e-2_R8, -1.e-5_R8, par_opt(6), istyle) ! enepar(ib,1) [ -0.01,  -0.0001]
      par_opt_phys(7) = fab(-1.e-1_R8, -1.e-3_R8, par_opt(7), istyle) ! enipar(ib,1) [ -0.10,  -0.001]
!     par_opt_phys(8) for separatrix shift
      par_opt_phys(8) = fab(-1.e-2_R8,  1.e-2_R8, par_opt(8), istyle)  ! shift [-0.01, 0.01]  [m]
      shift = par_opt_phys(8)

*   .. additional parameters for H-regime
      if(npar_opt.eq.16) then
       par_opt_phys(9)  =  fab(0.010_R8, 0.030_R8, par_opt(9), istyle)  !   delta1 [ 0.01,  0.03] [m]
       par_opt_phys(10) =  fab(0.001_R8, 0.005_R8, par_opt(10), istyle) !   delta2 [0.001, 0.005] [m]
       par_opt_phys(11) =  fab(0.050_R8, 0.300_R8, par_opt(11), istyle) !    alfa1  [0.05, 0.3] [m]
       par_opt_phys(12) =  fab(0.050_R8, 0.300_R8, par_opt(12), istyle) !    alfa2  [0.05, 0.3] [m]
       par_opt_phys(13) =  fab(1.000_R8, 3.000_R8, par_opt(13), istyle) !    alfa3  [1.00, 3.0] [m]
       par_opt_phys(14) =  fab(0.200_R8, 1.000_R8, par_opt(14), istyle) ! te_alfa1  [0.20, 1.0] [m]
       par_opt_phys(15) =  fab(0.200_R8, 1.000_R8, par_opt(15), istyle) ! te_alfa2  [0.20, 1.0] [m]
       par_opt_phys(16) =  fab(1.000_R8, 3.000_R8, par_opt(16), istyle) ! te_alfa3  [1.00, 3.0] [m]
      endif
      write(hlp_frm,'(a,i2,a)') '(a,',npar_opt,'es14.6)'
      write(71,hlp_frm) ' b2obfn: par_opt             ', par_opt
      write(71,hlp_frm) ' b2obfn: par_opt_phys        ', par_opt_phys
*   .. write binary file of physical parameters
      open(70,file='b2fparop.dat',form='unformatted')
      write(70) par_opt_phys
      close(70)

      call cpu_time(cpu_s)

#ifdef WINDOWS
      iresult = RUNQQ ('B2mnOpt','')
      inquire(file='b2fstate',exist=b2fstate_exist)
      if(b2fstate_exist) then
       lresult = SYSTEMQQ('del b2fstati')
       lresult = SYSTEMQQ('ren b2fstate b2fstati')
      endif
#else
      call system('./B2mnOpt.exe')
      inquire(file='b2fstate',exist=b2fstate_exist)
      if(b2fstate_exist) then
       call system('rm b2fstati')
       call system('cp b2fstate b2fstati')
      endif
#endif
      call cpu_time(cpu_f)
      write(71,"(' b2obfn takes ',es10.2,' sec')") cpu_f-cpu_s
* read the results of modelling
      call read_results
*   .. compute objective function and write the results of modelling
      func = 0.0_R8
      npoint=0
*   .. for inner lower
      call treat_exp_res(l_ne_inner_lower, fnm_ne_inner_lower,
     & ny, cy_ne_inner_lower, ne_inner_lower_mod, b_ne_inner_lower_mod,
     & c_ne_inner_lower_mod, d_ne_inner_lower_mod,
     & ndat_ne_inner_lower_exp, cy_ne_inner_lower_exp+shift,
     & ne_inner_lower_exp, npoint, func)

      call treat_exp_res(l_te_inner_lower, fnm_te_inner_lower,
     & ny, cy_te_inner_lower, te_inner_lower_mod,
     & b_te_inner_lower_mod, c_te_inner_lower_mod, d_te_inner_lower_mod,
     & ndat_te_inner_lower_exp, cy_te_inner_lower_exp+shift,
     & te_inner_lower_exp, npoint, func)

      call treat_exp_res(l_ti_inner_lower, fnm_ti_inner_lower,
     & ny, cy_ti_inner_lower, ti_inner_lower_mod,
     & b_ti_inner_lower_mod, c_ti_inner_lower_mod, d_ti_inner_lower_mod,
     & ndat_ti_inner_lower_exp, cy_ti_inner_lower_exp+shift,
     & ti_inner_lower_exp, npoint, func)

      call treat_exp_res(l_js_inner_lower, fnm_js_inner_lower,
     & ny, cy_js_inner_lower, js_inner_lower_mod,
     & b_js_inner_lower_mod, c_js_inner_lower_mod, d_js_inner_lower_mod,
     & ndat_js_inner_lower_exp,
     & cy_js_inner_lower_exp(1:ndat_js_inner_lower_exp)+shift,
     & js_inner_lower_exp, npoint, func)

*   .. for inner upper
      call treat_exp_res(l_ne_inner_upper, fnm_ne_inner_upper,
     & ny, cy_ne_inner_upper, ne_inner_upper_mod, b_ne_inner_upper_mod,
     & c_ne_inner_upper_mod, d_ne_inner_upper_mod,
     & ndat_ne_inner_upper_exp, cy_ne_inner_upper_exp+shift,
     & ne_inner_upper_exp, npoint, func)

      call treat_exp_res(l_te_inner_upper, fnm_te_inner_upper,
     & ny, cy_te_inner_upper, te_inner_upper_mod,
     & b_te_inner_upper_mod, c_te_inner_upper_mod, d_te_inner_upper_mod,
     & ndat_te_inner_upper_exp, cy_te_inner_upper_exp+shift,
     & te_inner_upper_exp, npoint, func)

      call treat_exp_res(l_ti_inner_upper, fnm_ti_inner_upper,
     & ny, cy_ti_inner_upper, ti_inner_upper_mod,
     & b_ti_inner_upper_mod, c_ti_inner_upper_mod, d_ti_inner_upper_mod,
     & ndat_ti_inner_upper_exp, cy_ti_inner_upper_exp+shift,
     & ti_inner_upper_exp, npoint, func)

      call treat_exp_res(l_js_inner_upper, fnm_js_inner_upper,
     & ny, cy_js_inner_upper, js_inner_upper_mod,
     & b_js_inner_upper_mod, c_js_inner_upper_mod, d_js_inner_upper_mod,
     & ndat_js_inner_upper_exp, cy_js_inner_upper_exp+shift,
     & js_inner_upper_exp, npoint, func)

*   .. for outer midplane
      call treat_exp_res(l_ne_outer_mdpl, fnm_ne_outer_mdpl,
     & ny, cy_ne_outer_mdpl, ne_outer_mdpl_mod,
     & b_ne_outer_mdpl_mod, c_ne_outer_mdpl_mod, d_ne_outer_mdpl_mod,
     & ndat_ne_outer_mdpl_exp, cy_ne_outer_mdpl_exp-shift,
     & ne_outer_mdpl_exp, npoint, func)

      call treat_exp_res(l_te_outer_mdpl, fnm_te_outer_mdpl,
     & ny, cy_te_outer_mdpl, te_outer_mdpl_mod,
     & b_te_outer_mdpl_mod, c_te_outer_mdpl_mod, d_te_outer_mdpl_mod,
     & ndat_te_outer_mdpl_exp, cy_te_outer_mdpl_exp-shift,
     & te_outer_mdpl_exp, npoint, func)

      call treat_exp_res(l_ti_outer_mdpl, fnm_ti_outer_mdpl,
     & ny, cy_ti_outer_mdpl, ti_outer_mdpl_mod,
     & b_ti_outer_mdpl_mod, c_ti_outer_mdpl_mod, d_ti_outer_mdpl_mod,
     & ndat_ti_outer_mdpl_exp, cy_ti_outer_mdpl_exp-shift,
     & ti_outer_mdpl_exp, npoint, func)

*   .. for outer lower
      call treat_exp_res(l_ne_outer_lower, fnm_ne_outer_lower,
     & ny, cy_ne_outer_lower, ne_outer_lower_mod, b_ne_outer_lower_mod,
     & c_ne_outer_lower_mod, d_ne_outer_lower_mod,
     & ndat_ne_outer_lower_exp, cy_ne_outer_lower_exp-shift,
     & ne_outer_lower_exp, npoint, func)

      call treat_exp_res(l_te_outer_lower, fnm_te_outer_lower,
     & ny, cy_te_outer_lower, te_outer_lower_mod,
     & b_te_outer_lower_mod, c_te_outer_lower_mod, d_te_outer_lower_mod,
     & ndat_te_outer_lower_exp, cy_te_outer_lower_exp-shift,
     & te_outer_lower_exp, npoint, func)

      call treat_exp_res(l_ti_outer_lower, fnm_ti_outer_lower,
     & ny, cy_ti_outer_lower, ti_outer_lower_mod,
     & b_ti_outer_lower_mod, c_ti_outer_lower_mod, d_ti_outer_lower_mod,
     & ndat_ti_outer_lower_exp, cy_ti_outer_lower_exp-shift,
     & ti_outer_lower_exp, npoint, func)

      call treat_exp_res(l_js_outer_lower, fnm_js_outer_lower,
     & ny, cy_js_outer_lower, js_outer_lower_mod,
     & b_js_outer_lower_mod, c_js_outer_lower_mod, d_js_outer_lower_mod,
     & ndat_js_outer_lower_exp, cy_js_outer_lower_exp-shift,
     & js_outer_lower_exp, npoint, func)

*   .. for outer upper
      call treat_exp_res(l_ne_outer_upper, fnm_ne_outer_upper,
     & ny, cy_ne_outer_upper, ne_outer_upper_mod, b_ne_outer_upper_mod,
     & c_ne_outer_upper_mod, d_ne_outer_upper_mod,
     & ndat_ne_outer_upper_exp, cy_ne_outer_upper_exp-shift,
     & ne_outer_upper_exp, npoint, func)

      call treat_exp_res(l_te_outer_upper, fnm_te_outer_upper,
     & ny, cy_te_outer_upper, te_outer_upper_mod,
     & b_te_outer_upper_mod, c_te_outer_upper_mod, d_te_outer_upper_mod,
     & ndat_te_outer_upper_exp, cy_te_outer_upper_exp-shift,
     & te_outer_upper_exp, npoint, func)

      call treat_exp_res(l_ti_outer_upper, fnm_ti_outer_upper,
     & ny, cy_ti_outer_upper, ti_outer_upper_mod,
     & b_ti_outer_upper_mod, c_ti_outer_upper_mod, d_ti_outer_upper_mod,
     & ndat_ti_outer_upper_exp, cy_ti_outer_upper_exp-shift,
     & ti_outer_upper_exp, npoint, func)

      call treat_exp_res(l_js_outer_upper, fnm_js_outer_upper,
     & ny, cy_js_outer_upper, js_outer_upper_mod,
     & b_js_outer_upper_mod, c_js_outer_upper_mod, d_js_outer_upper_mod,
     & ndat_js_outer_upper_exp, cy_js_outer_upper_exp-shift,
     & js_outer_upper_exp, npoint, func)

      if(l_power) then
       dif = (power_mod-power_exp)/power_exp
       func = func + dif**2
       npoint=npoint+1

       write(71,"(5x,'pw_mod',5x,'pw_exp',5x,'dif')")
       write(71,"(1x,3(es14.6))") power_mod, power_exp, dif
      endif

      func = sqrt(func/npoint)

      nfunc = nfunc + 1
      write(71,*) ' b2obfn: nfunc, func, npoint ', nfunc, func, npoint

* ..end of computation
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2obfn

      end subroutine b2obfn

      function fab(a, b, z, istyle)
      use b2mod_types
      implicit none
      integer istyle
      real (kind=R8) :: a, b, z, fab
      external subini, subend, xerrab

      call subini ('fab')
      select case (istyle)
        case(0)
!         x = z
          fab = z
        case(1)
!         x = b + (a-b)*exp(z)/(1+exp(z))
          if(z .le. 0.0_R8) then
            fab = b + (a-b)*exp(z)/(1.0_R8+exp(z))
          else
            fab = b + (a-b)/(1.0_R8+exp(-z))
          endif
        case(2)
!         x = (a+b)/2 + (b-a)/2*sin(z)
          fab = 0.5_R8*(a+b) + 0.5_R8*(b-a)*sin(z)
        case default
          write(71,*) ' fab: Incorrect istyle'
          call xerrab ('Incorrect istyle')
      end select
     
      call subend ()
      return
      end function fab

      function fcain(a, z, istyle)
      use b2mod_types
      implicit none
      integer istyle
      real (kind=R8) :: a, z, fcain
      external subini, subend, xerrab

      call subini ('fcain')
      select case (istyle)
        case(0)
!         x = z
          fcain = z
        case(1)
!         x = a + exp(z)
          fcain = a + exp(z)
        case(2)
!         x = a + z**2
          fcain = a + z**2
        case default
          write(71,*) ' fcain: Incorrect istyle'
          call xerrab ('Incorrect istyle')
      end select
      call subend ()
      return
      end function fcain

      subroutine treat_exp_res(l_ne_outer_mdpl, fnm_ne_outer_mdpl,
     &  ny, cy_ne_outer_mdpl, ne_outer_mdpl_mod,
     &  b_ne_outer_mdpl_mod, c_ne_outer_mdpl_mod, d_ne_outer_mdpl_mod,
     &  ndat_ne, cy_ne_outer_mdpl_exp, ne_outer_mdpl_exp, npoint, func)
      use b2mod_types
      implicit none
      logical(1) :: l_ne_outer_mdpl
      character(*) :: fnm_ne_outer_mdpl
      integer :: ny, npoint, ndat_ne
      real (kind=R8) :: func,
     &   cy_ne_outer_mdpl(-1:ny),
     &   ne_outer_mdpl_mod(-1:ny), b_ne_outer_mdpl_mod(-1:ny),
     & c_ne_outer_mdpl_mod(-1:ny), d_ne_outer_mdpl_mod(-1:ny),
     & cy_ne_outer_mdpl_exp(1:ndat_ne), ne_outer_mdpl_exp(1:ndat_ne)

      integer :: iy
      real (kind=R8) :: dif, nemod, dseval
      character(128) :: fmt

      fmt="(8x,'y',7x,1x,'"//trim(fnm_ne_outer_mdpl)//"_exp'"//
     &  ",4x,'"//trim(fnm_ne_outer_mdpl)//"_mod',8x,'dif')"
      if(l_ne_outer_mdpl) then
        call approx_profile_mod(ny+2,
     &                 cy_ne_outer_mdpl, ne_outer_mdpl_mod,
     &   b_ne_outer_mdpl_mod, c_ne_outer_mdpl_mod, d_ne_outer_mdpl_mod)
        write(71,trim(fmt))
        do iy=1,ndat_ne
          if(cy_ne_outer_mdpl(-1)*1.0001.le.cy_ne_outer_mdpl_exp(iy)
     &        .and.
     &     cy_ne_outer_mdpl_exp(iy).le.cy_ne_outer_mdpl(ny)*1.0001) then
            nemod = dseval(ny+2, cy_ne_outer_mdpl_exp(iy),
     &         cy_ne_outer_mdpl, ne_outer_mdpl_mod, b_ne_outer_mdpl_mod,
     *         c_ne_outer_mdpl_mod, d_ne_outer_mdpl_mod)
            dif = (nemod-ne_outer_mdpl_exp(iy))/ne_outer_mdpl_exp(iy)
            func = func + dif**2
            npoint=npoint+1
            write(71,"(1x,4(es14.6))") cy_ne_outer_mdpl_exp(iy),
     &              ne_outer_mdpl_exp(iy), nemod, dif
          endif
        enddo
      endif
      end subroutine treat_exp_res

