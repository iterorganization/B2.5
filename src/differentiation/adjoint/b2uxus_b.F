*-----------------------------------------------------------------------
*.specification

      subroutine b2uxus_b(nCv, mpg, aa, aab, itcnt,   
     .  resfun, resfunb, corfun, corfunb, name)
      use b2mod_types
      use b2mod_ma28_for_us
      use b2us_map_diff
      use b2mod_lwmain
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, itcnt
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   aa(13*nCv), resfun(nCv)
      real (kind=R8) ::
     *   aab(13*nCv), resfunb(nCv)
      character*(*) name 
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   corfun(nCv)
      real (kind=R8) ::
     *   corfunb(nCv)
*   ..common blocks

*-----------------------------------------------------------------------
*     Rearranged routine which solves a linear system of equation A*x=b
*     in adjoint mode.
*     AT*bb=xb --> bb = (A^-1)^T*xb
*     Ab(i,j) = -x(j)*bb(i)
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nCv - integer, input.
*     See the calling routine for description.
*
*     aa - () real*8 array, input.
*     aa contains the stencil.
*
*     itcnt - integer, input.
*     itcnt specifies the iteration counter in the calling routine.
*
*     resfun - (nCv) real*8 array, input.
*     resfun specifies the right hand side.
*
*     corfun - (nCv) real*8 array, output.
*     corfun contains the computed solution.
*
*     (parameters in common:)
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, style, itmax, iters, itpar, cpu
      parameter (itmax=100)
      real (kind=R8) ::
     *   acpar, nrmc(0:itmax-1), cpu_start, cpu_end, c_end, c_start
      integer, save :: nz
      integer, allocatable, save :: ind(:), iasave(:), jasave(:)
      real (kind=R8) :: incrbb(nCv)
      integer, allocatable :: dummy(:)
      integer :: iCv, iNv
      integer mult_nonzero, mult_solvdim, mult_solvdim1
*   ..procedures
      intrinsic min, cpu_time
      external xertst, ipgeti, ipgetr,
     &  slv5pt, iluter, xerrab
*   ..initialization
      save ncall, style, acpar, cpu
      data ncall/0/, style/2/, acpar/8.0e0_R8/, cpu/0/
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2uxus_b')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti ('b2ux9p_style', style)
       call xertst (0.le.style.and.style.lt.3,
     &   'faulty internal parameter style')
       call ipgetr ('b2ux9p_acpar', acpar)
       call xertst (0.lt.acpar, 'faulty internal parameter acpar')
       call ipgeti ('b2ux9p_cpu', cpu)
       if(style.eq.2) then
         mult_nonzero=10
         mult_solvdim=15
         mult_solvdim1=10
         call ipgeti ('b2uxus_mult_nonzero', mult_nonzero)
         call ipgeti ('b2uxus_mult_solvdim', mult_solvdim)
         call ipgeti ('b2uxus_mult_solvdim1', mult_solvdim1)
         nnd=nCv
         n_nonzeros=mpg%nCmxNv
         solvdim=mult_solvdim*n_nonzeros
         solvdim1=mult_solvdim1*n_nonzeros
         call alloc_b2mod_ma28_for_us                                  !sv 13.10.06
         irnh = 0                                                      !sv 15.05.07
         icnh = 0                                                      !sv 15.05.07
         nltrsol=2
         call ipgeti ('b2uxus_nltrsol', nltrsol)
       endif
      endif
*   ..test nCv, itcnt
      call xertst (0.le.nCv, 'faulty argument nCv')
      call xertst (0.le.itcnt, 'faulty argument itcnt')
*   ..obtain work space
      if (ncall.eq.0) then
       write(*,*) 'Allocating space for solver b2uxus_b'
       if (style.eq.0) then
        call xerrab ('b2uxus_b: Style 0 is not implemented')
      else if (style.eq.1) then
        call xerrab ('b2uxus_b: Style 1 is not implemented')
      else if (style.eq.2) then
        allocate(ind(1:solvdim))
        allocate(iasave(1:solvdim),jasave(1:solvdim))
        allocate(dummy(1:solvdim))
      endif
      write(*,*) 'Allocated space for solver b2uxus'
      endif

*   ..solve the correction equation
      call cpu_time(cpu_start)
      if (style.eq.0) then
*    ..use iluter
       call xerrab ('b2uxus: Style 0 is not implemented')
      else if (style.eq.1) then
*    ..use slv5pt
       call xerrab ('b2uxus: Style 1 is not implemented')
      else if (style.eq.2) then
       call cpu_time(c_start)
       call storage_us (nCv, mpg, aa, n_nonzeros, irnh, icnh, a1h, nz,
     *                 ind, iasave, jasave, name)       
       call cpu_time(c_end)
*csc  Swap indexes to have transpose matrix
       dummy = irnh
       irnh = icnh
       icnh = dummy 
       if(cpu.gt.0) write(*,'(a,f6.3)') 'cpu for storage_us ',
     &               c_end-c_start
*    ..use ma28
       call cpu_time(cpu_start)

*csc  Solve the system AT*bb=xb / AT*incrbb=corfunb
        call ILUTERN_US(NCV, NZ, corfunb, incrbb, name)
        resfunb = resfunb + incrbb
      endif

*csc  Obtain x, solution of the original system A*x=b
      call b2uxus (nCv, mpg, aa, itcnt, resfun, corfun, name)
*     Obtain Ab(i,j) = -x(j)*bb(i)
      do iCv=1,nCv
        do iNv=1,mpg%cvNvP(iCv,2)
          aab(mpg%cvNvP(iCv,1)+iNv-1) = aab(mpg%cvNvP(iCv,1)+iNv-1) - 
     &      corfun(mpg%cvNv(mpg%cvNvP(iCv,1)+iNv-1))*incrbb(iCv)
        enddo
        corfunb(iCv) = 0.0
      enddo

      call cpu_time(cpu_end)
      if(cpu.gt.0) write(*,'(a,f6.3)') 'cpu for '//trim(name),
     &             cpu_end-cpu_start

* ..return
      ncall = ncall+1
      call subend ()
      return
         
*-----------------------------------------------------------------------
*.end b2uxus_b

      end subroutine b2uxus_b



