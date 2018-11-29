*-----------------------------------------------------------------------
*.specification

      subroutine b2uxus (nCv, mpg, aa, itcnt, resfun, corfun, name)    !sv 13.10.06 {
      use b2mod_types
      use b2mod_ma28_for_us
      use b2us_map
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, itcnt
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   aa(13*nCv), resfun(nCv)
      character*(*) name                                              !sv 18.05.02
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   corfun(nCv)
*   ..common blocks

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2UXUS solves a system of equations on an unstructured mesh. 
*     This routine is an interface to an iterative or a direct solver;
*     the choice depends upon an internal parameter.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine solves a system of equations of the form associated
*     with the five-point discretization of a single scalar partial
*     differential equation on a topologically rectangular mesh. The
*     routine is an interface to one of two solvers. If the internal
*     parameter "style" has the value 0 then a solver based on Stone's
*     SIP-ILU procedure is used; if style.eq.1 then a direct band
*     solver is used.
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
      real (kind=R8), allocatable, save ::
     1 aaf(:),wrk0(:),wrk1(:),wrk2(:)
      integer, allocatable, save :: iwrk(:)

!      real (kind=R8), allocatable :: a(:)
      integer, save :: nz
!      integer, allocatable :: ia(:), ja(:)
      integer, allocatable, save :: ind(:), iasave(:), jasave(:)

      integer mult_nonzero, mult_solvdim, mult_solvdim1
*   ..procedures
      intrinsic min
      external subini, subend, xertst, ipgeti, ipgetr, lwget, lwiget,
     &  slv5pt, iluter, xerrab
*   ..initialization
      save ncall, style, acpar, cpu
      data ncall/0/, style/2/, acpar/8.0e0_R8/, cpu/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2uxus')
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
       write(*,*) 'Allocating space for solver b2uxus'
       if (style.eq.0) then
        call xerrab ('b2uxus: Style 0 is not implemented')
      else if (style.eq.1) then
        call xerrab ('b2uxus: Style 1 is not implemented')
      else if (style.eq.2) then
        allocate(ind(1:solvdim))
        allocate(iasave(1:solvdim),jasave(1:solvdim))
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
       if(cpu.gt.0) write(*,'(a,f6.3)') 'cpu for storage_us ',
     &               c_end-c_start
*    ..use ma28
       call cpu_time(cpu_start)
       call ILUTERN_US(NCV, NZ, resfun, corfun, name)     !sv 13.10.06
      endif
      call cpu_time(cpu_end)
      if(cpu.gt.0) write(*,'(a,f6.3)') 'cpu for '//trim(name),
     &             cpu_end-cpu_start

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2uxus

      end

      subroutine storage_us (nCv, mpg, aa, nnz, ia, ja, a, k,            !sv 22.11.06
     *                      ind, iasave, jasave, name)
      use b2mod_types
      use b2us_map
      implicit none
*  ..input
      integer nCv, nnz
      type (mapping), intent (in) :: mpg
      real (kind=R8) :: aa(13*nCv)
      character*(*) name
*  ..output
      integer ia(nnz), ja(nnz), k, ind(nnz), iasave(nnz), jasave(nnz)
      real (kind=R8) :: a(nnz)
!                 2_____3_____4                     5___________6
!                 |     |     |                     |     |     |
!                 |     |     |                     |     |     |
!                -1_____0_____1                     |_____|_____|
!                 |     |     |                     |     |     |
!                 |     |     |                     |     |     |
!                -4____-3____-2                    -6_____|____-5

!                (ix-1,iy+1)______(ix,iy+1)____(ix+1,iy+1)
!                     |              |              |
!                     |              |              |
!                     |              |              |
!                (ix-1,iy  )______(ix,iy)______(ix+1,iy  )
!                     |              |              |
!                     |              |              |
!                     |              |              |
!                (ix-1,iy-1)______(ix,iy-1)____(ix+1,iy-1)
*  ..locals
      logical :: lord
      logical, save :: lfirst=.true.
      integer :: i, iCv, iNv, ner, jtemp, ipt
      integer :: ip(1:nCv+1)
      real (kind=R8) :: temp
      integer, save :: ncall = 0
      integer, save :: kwrite_fmt = -1
      integer, save :: kwrite_bnr = -1

      if(ncall.eq.0) then
        call ipgeti ('storage_us_kwrite_fmt', kwrite_fmt)
        call ipgeti ('storage_us_kwrite_bnr', kwrite_bnr)
      endif

      if (lfirst) then
        ind=0
        k=0
        ip(1)=1
        do iCv = 1,nCv
          ner=0
          do iNv = 1,mpg%cvNvP(iCv,2)
            if (aa(mpg%cvNvP(iCv,1)+iNv-1) .ne. 0.0e0_R8 ) then
              k      = k+1
              ner    = ner+1
              a(k)   = aa(mpg%cvNvP(iCv,1)+iNv-1)
              ind(k) = mpg%cvNvP(iCv,1)+iNv-1
              ia(k)  = iCv
              ja(k)  = mpg%cvNv(mpg%cvNvP(iCv,1)+iNv-1)
            endif
          enddo
          ip(iCv+1)=ip(iCv) + ner
          do
            lord=.true.
            do ipt=ip(iCv),ip(iCv+1)-2
              if(ja(ipt).gt.ja(ipt+1)) then
                jtemp     = ja(ipt)
                ja(ipt)   = ja(ipt+1)
                ja(ipt+1) = jtemp
                temp     = a(ipt)
                a(ipt)   = a(ipt+1)
                a(ipt+1) = temp
                jtemp      = ind(ipt)
                ind(ipt)   = ind(ipt+1)
                ind(ipt+1) = jtemp
                lord=.false.
              endif
            enddo
            if (lord) exit
          enddo
        enddo

        do  i=1,k
          iasave(i)=ia(i)
          jasave(i)=ja(i)
        enddo                                              !sv 15.05.07 }
!        lfirst=.false.

      else
        do i=1,k
          ia(i)=iasave(i)                                   !sv 15.05.07 }
          ja(i)=jasave(i)
          a(i)=aa(ind(i))
        enddo
      endif

      if(kwrite_bnr.eq.ncall) then
!        write(*,*) 'Opening 76'
!        open(76,file=
!     &   'st_9_rhs_binary_'//trim(name)//'.dat',form='unformatted')
!         write(76) n
!         do i=1,n
!          write(76) r1(i)
!         enddo
!        close(76)
        write(*,*) 'Opening 77'
        open(77,file=
     &   'st_us_matrix_binary_'//trim(name)//'.dat',form='unformatted') !srv 11.03.13
        write(77) nCv, k
        do i=1,k
          write(77) ia(i), ja(i), a(i)
        enddo
        close(77)
      endif
      if(kwrite_fmt.eq.ncall) then
!        write(*,*) 'Opening 76'
!        open(76,file='st_9_rhs_'//trim(name)//'.dat')
!         write(76,*) n
!         do i=1,n
!          write(76,*) r1(i)
!         enddo
!        close(76)
        open(77,file='st_us_matrix_'//trim(name)//'.dat')
        write(77,*) nCv, k
        do i=1,k
          write(77,*) ia(i), ja(i), a(i)
        enddo
        close(77)
      endif

      ncall = ncall + 1

      end                                                             !sv 22.11.06

