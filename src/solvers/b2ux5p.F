*-----------------------------------------------------------------------
*.specification

      subroutine b2ux5p (nx, ny, aa5, itcnt, resfun, corfun, name)   !srv 18.03.10
      use b2mod_types
      use b2mod_ma28
      use b2mod_indirect
      use b2mod_openmp
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, itcnt
      character*(*) name                                             !srv 18.03.10
      real (kind=R8) ::
     *   aa5(-1:nx,-1:ny,-2:2), resfun(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   corfun(-1:nx,-1:ny)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2UX5P solves a five-point system of equations. This routine is
*     an interface to an iterative or a direct solver; the choice
*     depends upon an internal parameter.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine solves a system of equations of the form associated
*     with the five-point discretisation of a single scalar partial
*     differential equation on a topologically rectangular mesh. The
*     routine is an interface to one of two solvers. If the internal
*     parameter "style" has the value 0 then a solver based on the Stone
*     SIP-ILU procedure is used; if style.eq.1 then a direct band
*     solver is used.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nx, ny - integer, input.
*     See the calling routine for description.
*
*     aa5 - (-1:nx,-1:ny,-2:2) real*8 array, input.
*     aa5 contains the five-point matrix.
*
*     itcnt - integer, input.
*     itcnt specifies the iteration counter in the calling routine.
*
*     resfun - (-1:nx,-1:ny) real*8 array, input.
*     resfun specifies the right hand side.
*
*     corfun - (-1:nx,-1:ny) real*8 array, output.
*     corfun contains the computed solution.
*
*     (parameters in common:)
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, style, itmax, iters, itpar, cpu
      parameter (itmax=100)
#ifdef _OPENMP
      real (kind=R8) :: cpu_start, cpu_end
#else
      real (kind=R4) :: cpu_start, cpu_end
#endif
      real (kind=R8) :: acpar, nrmc(0:itmax-1)
*   ..procedures
      intrinsic min
      external xertst, ipgeti, ipgetr, slv5pt, iluter
#if ( defined(_OPENMP) && defined(CRAY) )
      real (kind=R8) :: omp_get_wtime
      external omp_get_wtime
#endif
*   ..initialisation
      save ncall, style, acpar, cpu
#ifdef _OPENMP
!!FR The save variables used by ilutern need to be declared threadprivate
!!FR to ensure that each thread accesses its own copy of these and not a
!!FR global version. The iluter and slv5pt specific arrays have also been
!!FR declared threadprivate to ensure any attempts to allocate them on
!!FR multiple threads will work. However, the code will be stopped if
!!FR the iluter or slv5pt are used as these have not been parallelised.
!$OMP THREADPRIVATE (ncall, style, acpar, cpu)
#endif
      data ncall/0/, style/2/, acpar/8.0_R8/, cpu/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2ux5p')
*   ..test nx, ny, itcnt
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (0.le.itcnt, 'faulty argument itcnt')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti ('b2ux5p_style', style)
       call xertst (0.le.style.and.style.lt.3,
     &   'faulty internal parameter style')
       call ipgetr ('b2ux5p_acpar', acpar)
       call xertst (0.lt.acpar, 'faulty internal parameter acpar')
       call ipgeti ('b2ux5p_cpu', cpu)
      endif

*   ..solve the correction equation
#ifdef _OPENMP
      cpu_start = omp_get_wtime()
#else
      call cpu_time(cpu_start)
#endif
      if (style.eq.0) then
*    ..use iluter
       iters = min(itmax,4)
       itpar = itcnt*iters
!!FR iluter and slv5pt are not parallelised thus we need to stop the
!!FR code if attempting to run them in parallel
       if (in_parallel()) then
        stop 'STOP: iluter not parallelised, re-run in serial'
       else
        call iluter (nx+2, ny+2, aa5, resfun, iters, acpar, itpar,
     &   corfun, nrmc, aaf_5p, wrk0_5p, wrk1_5p, wrk2_5p)
       endif
      else if (style.eq.1) then
*    ..use slv5pt
       iters = 1
       if (in_parallel()) then
        stop 'STOP: slv5pt not parallelised, re-run in serial'
       else
        call slv5pt (nx+2, ny+2, aa5, resfun, iters, corfun, nrmc,
     &   min(nx+2,ny+2), aaf_5p, wrk0_5p, iwrk_5p)
       endif
      else if (style.eq.2) then
*    ..use ma28
        call ilutern(NX+2, NY+2, NX+2,
     1   aa5(-1,-1,-2), aa5(-1,-1,-1), aa5(-1,-1,0),
     2   aa5(-1,-1,1), aa5(-1,-1,2),
     .   resfun, corfun, name)                   !srv 18.03.10
      endif
#ifdef _OPENMP
      cpu_end = omp_get_wtime()
#else
      call cpu_time(cpu_end)
#endif
      if(cpu.gt.0) write(*,'(a,f6.3)') 'cpu ',cpu_end-cpu_start

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2ux5p

      end subroutine b2ux5p
