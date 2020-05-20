*-----------------------------------------------------------------------
*.specification

      subroutine b2ux5p (nx, ny, aa5, itcnt, resfun, corfun, name)   !srv 18.03.10
      use b2mod_types
      use b2mod_ma28
      use b2mod_indirect
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, itcnt
      character*(*) name                                             !srv 18.03.10
      real (kind=R8) ::
     *   aa5(-1:nx,-1:ny,-2:2), resfun(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   corfun(-1:nx,-1:ny)
*   ..common blocks

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
      integer ncall, style, n2, itmax, iters, itpar, cpu, ix, iy
      parameter (itmax=100)
      real (kind=R8) ::
     *   acpar, nrmc(0:itmax-1), cpu_start, cpu_end
      real (kind=R8), allocatable, save ::
     1 aaf(:),wrk0(:),wrk1(:),wrk2(:)
      integer, allocatable, save :: iwrk(:)
      integer mult_nonzero, mult_solvdim, mult_solvdim1
*   ..procedures
      intrinsic min, cpu_time
      external subini, subend, xertst, ipgeti, ipgetr, slv5pt, iluter
*   ..initialisation
      save ncall, style, acpar, cpu
      data ncall/0/, style/2/, acpar/8.0_R8/, cpu/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2ux5p')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti ('b2ux5p_style', style)
       call xertst (0.le.style.and.style.lt.3,
     &   'faulty internal parameter style')
       call ipgetr ('b2ux5p_acpar', acpar)
       call xertst (0.lt.acpar, 'faulty internal parameter acpar')
       call ipgeti ('b2ux5p_cpu', cpu)
       if(style.eq.2) then
         mult_nonzero=10
         mult_solvdim=15
         mult_solvdim1=10
         call ipgeti ('b2ux5p_mult_nonzero', mult_nonzero)
         call ipgeti ('b2ux5p_mult_solvdim', mult_solvdim)
         call ipgeti ('b2ux5p_mult_solvdim1', mult_solvdim1)
         nnd=(nx+2)*(ny+2)
         n_nonzeros=mult_nonzero*nnd
         solvdim=mult_solvdim*n_nonzeros
         solvdim1=mult_solvdim1*n_nonzeros
         call alloc_b2mod_ma28
         call iniilu(nx+2,ny+2)
         nltrsol=2
         call ipgeti ('b2ux5p_nltrsol', nltrsol)
       endif
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, itcnt
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (0.le.itcnt, 'faulty argument itcnt')
*   ..obtain work space
      if (ncall.eq.0) then
        write(*,*) 'Allocating space for solver'
        if (style.eq.0) then
*    ..obtain workspace for iluter
          allocate(aaf(1:n2*5))
          allocate(wrk0(1:n2))
          allocate(wrk1(1:n2))
          allocate(wrk2(1:n2))
          allocate(iwrk(1:1))
        else if (style.eq.1) then
*    ..obtain workspace for slv5pt
          allocate(aaf(1:n2*(3*min(nx+2,ny+2)+1)))
          allocate(wrk0(1:n2))
          allocate(wrk1(1:1))
          allocate(wrk2(1:1))
          allocate(iwrk(1:n2))
        else if (style.eq.2) then
        endif
        if (style.lt.2) then
          do ix = -1, nx
            do iy = -1, ny
              call xertst (
     .    (leftix(ix,iy).eq.-2 .or. leftix(ix,iy).eq.ix-1) .and.
     .     leftiy(ix,iy).eq.iy .and.
     .    (rightix(ix,iy).eq.nx+1 .or. rightix(ix,iy).eq.ix+1) .and.
     .     rightiy(ix,iy).eq.iy .and.
     .     bottomix(ix,iy).eq.ix .and.
     .    (bottomiy(ix,iy).eq.-2 .or. bottomiy(ix,iy).eq.iy-1) .and.
     .     topix(ix,iy).eq.ix .and.
     .    (topiy(ix,iy).eq.ny+1 .or. topiy(ix,iy).eq.iy+1),
     .   'Solver type inappropriate for non-trivial neighbourhoods')
            enddo
          enddo
        endif
        write(*,*) 'Allocated space for solver'
      endif

*   ..solve the correction equation
      call cpu_time(cpu_start)
      if (style.eq.0) then
*    ..use iluter
       iters = min(itmax,4)
       itpar = itcnt*iters
       call iluter (nx+2, ny+2, aa5, resfun, iters, acpar, itpar,
     &   corfun, nrmc, aaf, wrk0, wrk1, wrk2)
      else if (style.eq.1) then
*    ..use slv5pt
       iters = 1
       call slv5pt (nx+2, ny+2, aa5, resfun, iters, corfun, nrmc,
     &   min(nx+2,ny+2), aaf, wrk0, iwrk)
      else if (style.eq.2) then
*    ..use ma28
        call ilutern(NX+2, NY+2, NX+2,
     1   aa5(-1,-1,-2), aa5(-1,-1,-1), aa5(-1,-1,0),
     2   aa5(-1,-1,1), aa5(-1,-1,2),
     .   resfun, corfun, name)                   !srv 18.03.10
      endif
      call cpu_time(cpu_end)
      if(cpu.gt.0) write(*,'(a,f6.3)') 'cpu ',cpu_end-cpu_start

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2ux5p

      end subroutine b2ux5p
