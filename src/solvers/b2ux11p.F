*-----------------------------------------------------------------------
*.specification

      subroutine b2ux11p (nx, ny, aa11, itcnt, resfun, corfun, name)    !srv 02.04.18 {
      use b2mod_types
      use b2mod_illtern_7
      use b2mod_ilutern_11
      use b2mod_ma28_for_11diag
      use b2mod_ysmp_sdrv
      use b2mod_indirect
      use b2mod_openmp
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, itcnt
      real (kind=R8) ::
     *   aa11(-1:nx,-1:ny,-7:7), resfun(-1:nx,-1:ny)
      character*(*) name                                              !srv 18.05.02
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   corfun(-1:nx,-1:ny)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2UX11P solves an eleven-point system of equations. This routine is
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
*     parameter "style" has the value 0 then a solver based on Stone
*     SIP-ILU procedure is used; if style.eq.1 then a direct band
*     solver is used.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nx, ny - integer, input.
*     See the calling routine for description.
*
*     aa11 - (-1:nx,-1:ny,-7:7) real*8 array, input.
*     aa11 contains x,y-part of the 11-point stencil.
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
      integer ncall, style, n2, itmax, iters, itpar, cpu
      integer, save :: nz
      parameter (itmax=100)
#ifdef _OPENMP
      real (kind=R8) :: cpu_start, cpu_end, c_end, c_start
#else
      real (kind=R4) :: cpu_start, cpu_end, c_end, c_start
#endif
      real (kind=R8) :: acpar, nrmc(0:itmax-1)
      real (kind=R8), allocatable, save ::
     1 aaf(:),wrk0(:),wrk1(:),wrk2(:)
      integer, allocatable, save :: iwrk(:)
      integer, allocatable, save :: indices(:,:), ipsave(:), jpsave(:)
      integer mult_nonzero, mult_solvdim, mult_solvdim1
*   ..procedures
      intrinsic min
      external xertst, ipgeti, ipgetr, slv5pt, iluter, xerrab
#if ( defined(_OPENMP) && defined(CRAY) )
      real (kind=R8) :: omp_get_wtime
      external omp_get_wtime
#endif
*   ..initialisation
      save ncall, style, acpar, cpu
      !$OMP THREADPRIVATE (iwrk, aaf, wrk0, wrk1, wrk2)
      !$OMP THREADPRIVATE (ncall, style, acpar, cpu, nz)
      !$OMP THREADPRIVATE (indices, ipsave, jpsave)
      data ncall/0/, style/2/, acpar/8.0e0_R8/, cpu/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2ux11p')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti ('b2ux11p_style', style)
       call xertst (0.le.style.and.style.lt.4,
     &   'faulty internal parameter style')
       call ipgetr ('b2ux11p_acpar', acpar)
       call xertst (0.lt.acpar, 'faulty internal parameter acpar')
       call ipgeti ('b2ux11p_cpu', cpu)
       if(style.eq.2) then
         mult_nonzero=12
         mult_solvdim=15
         mult_solvdim1=10
         call ipgeti ('b2ux11p_mult_nonzero', mult_nonzero)
         call ipgeti ('b2ux11p_mult_solvdim', mult_solvdim)
         call ipgeti ('b2ux11p_mult_solvdim1', mult_solvdim1)
         nnd=(nx+2)*(ny+2)
         n_nonzeros=mult_nonzero*nnd
         solvdim=mult_solvdim*n_nonzeros
         solvdim1=mult_solvdim1*n_nonzeros
         call alloc_b2mod_ma28_for_11diag                          !srv 02.04.18
!         call iniilu_9(nx+2,ny+2)                                 !srv 13.10.06
         irnh = 0                                                  !srv 15.05.07
         icnh = 0                                                  !srv 15.05.07
         nltrsol=2
         call ipgeti ('b2ux11p_nltrsol', nltrsol)
       endif
       if(style.eq.3) then                                         !srv 29.05.02 {
         nnd=(nx+2)*(ny+2)
         nnd_sdrv  = (nx+2)*(ny+2)
         n_nz_sdrv = 4*nnd_sdrv
!         nsp = 3*nnd + 12*n_nz_sdrv  ! was 7*n_nz_sdrv without vis per !srv 13.02.03
         nsp = 3*nnd + 13*n_nz_sdrv                                !srv 20.02.06
         call alloc_b2mod_ysmp_sdrv
         ia = 0                                                    !srv 15.05.07
         ja = 0                                                    !srv 15.05.07
       endif                                                       !srv 29.05.02 }
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, itcnt
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (0.le.itcnt, 'faulty argument itcnt')
*   ..obtain work space
      if (ncall.eq.0) then
        write(*,*)
     .   'Allocating space for solver for 11-point-stencil matrix'  !srv 11.01.00
        if (style.eq.0) then
          call xerrab ('Style not applicable to 11-point-stencil!')
*    ..obtain workspace for iluter
          allocate(aaf(1:n2*5))
          allocate(wrk0(1:n2))
          allocate(wrk1(1:n2))
          allocate(wrk2(1:n2))
          allocate(iwrk(1:1))
        else if (style.eq.1) then
          call xerrab ('Style not applicable to 11-point-stencil!')
*    ..obtain workspace for slv5pt
          allocate(aaf(1:n2*(3*min(nx+2,ny+2)+1)))
          allocate(wrk0(1:n2))
          allocate(wrk1(1:1))
          allocate(wrk2(1:1))
          allocate(iwrk(1:n2))
        else if (style.eq.2) then
          allocate(indices(1:3,1:solvdim))
          allocate(ipsave(1:solvdim),jpsave(1:solvdim))
        else if (style.eq.3) then
          allocate(ipsave(1:n_nz_sdrv))
        endif
        write(*,*)
     .   'Allocated space for solver for 11-point-stencil matrix'   !srv 11.01.00
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
       call xerrab ('b2ux11p: Style 0 is not implemented')
       if (in_parallel()) then
        stop 'STOP: iluter not parallelised, re-run in serial'
       else
        call iluter (nx+2, ny+2, aa11, resfun, iters, acpar, itpar,
     &   corfun, nrmc, aaf, wrk0, wrk1, wrk2)
       endif
      else if (style.eq.1) then
*    ..use slv5pt
       iters = 1
       call xerrab ('b2ux11p: Style 1 is not implemented')
       if (in_parallel()) then
        stop 'STOP: slv5pt not parallelised, re-run in serial'
       else
        call slv5pt (nx+2, ny+2, aa11, resfun, iters, corfun, nrmc,
     &   min(nx+2,ny+2), aaf, wrk0, iwrk)
       endif
      else if (style.eq.2) then
#ifdef _OPENMP
       c_start = omp_get_wtime()
#else
       call cpu_time(c_start)
#endif
       call storage_11 (nx, ny, aa11, n_nonzeros, irnh, icnh, a1h, nz,
     *                  indices, ipsave, jpsave, name)
#ifdef _OPENMP
       c_end = omp_get_wtime()
#else
       call cpu_time(c_end)
#endif
       if(cpu.gt.0) write(*,'(a,f6.3)') 'cpu for storage_11 ',
     &               c_end-c_start
*    ..use ma28
#ifdef _OPENMP
       cpu_start = omp_get_wtime()
#else
       call cpu_time(cpu_start)
#endif
!       call ilutern_9(NX+2, NY+2, NX+2,
!    &   aa11(-1,-1,-3),
!    &   aa11(-1,-1,-2), aa11(-1,-1,-1), aa11(-1,-1,0),
!    &   aa11(-1,-1,1), aa11(-1,-1,2), aa11(-1,-1,3),
!    .   resfun, corfun, name)
       call ILUTERN_11(NX+2, NY+2, NX+2, NZ, resfun, corfun, name)  !srv 02.04.18
      else if (style.eq.3) then                                     !srv 29.05.02 {
!        call cpu_time(cpu_start)
        call storage_11 (nx, ny, aa11, n_nz_sdrv, ia, ja, a, nz,
     *                   indexes, ipsave, jasave, name)
!        call cpu_time(cpu_end)
!        write (*,"(' Time of transform is',1pd9.2)") cpu_end-cpu_start
!        write(*,*) 'STORAGE_11: nz = ',nz
        call xertst (nz.le.N_NZ_SDRV, ' nz.gt.N_NZ_SDRV')
*    ..use sdrv from YSMP
!        call cpu_time(cpu_start)
        call illtern_7(NX+2, NY+2, NX+2, NZ, resfun, corfun)        !srv 29.05.02 }
!        call cpu_time(cpu_end)
!        write (*,"(' Time of illtern_7 is',1pd9.2)") cpu_end-cpu_start
      endif
#ifdef _OPENMP
      cpu_end = omp_get_wtime()
#else
      call cpu_time(cpu_end)
#endif
      if(cpu.gt.0) write(*,'(a,f6.3)') 'cpu for '//trim(name),
     &                                  cpu_end-cpu_start

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2ux11p

      end subroutine b2ux11p                                           !srv 02.04.18 }

      subroutine storage_11 (nx, ny, aa11, nnz, ia, ja, a, k,          !srv 02.04.18 {
     *                       indexes, iasave, jasave, name)
      use b2mod_types
      use b2mod_openmp
      use b2mod_indirect
      implicit none
*  ..input
      integer nx, ny, nnz
      real (kind=R8) :: aa11(-1:nx,-1:ny,-7:7)
      character*(*) name
*  ..output
      integer ia(nnz), ja(nnz), k,
     *        indexes(1:3,nnz), iasave(nnz), jasave(nnz)
      real (kind=R8) :: a(nnz)
!                                                           7
!                                                           |
!                                                           |
!       2_____3_____4          5___________6           _____3_____
!       |     |     |          |     |     |          |     |     |
!       |     |     |          |     |     |          |     |     |
!      -1_____0_____1          |_____|_____|          |_____0_____|
!       |     |     |          |     |     |          |     |     |
!       |     |     |          |     |     |          |     |     |
!      -4____-3____-2         -6_____|____-5          |____-3_____|
!                                                           |
!                                                           |
!                                                          -7
!
!                                 (ix,iy+2)
!                                    |
!                                    |
!                                    |
!                (ix-1,iy+1)______(ix,iy+1)____(ix+1,iy+1)
!                     |              |              |
!                     |              |              |
!                     |              |              |
!                (ix-1,iy  )______(ix,iy)______(ix+1,iy  )
!                     |              |              |
!                     |              |              |
!                     |              |              |
!                (ix-1,iy-1)______(ix,iy-1)____(ix+1,iy-1)
!                                    |
!                                    |
!                                    |
!                                (ix,iy-2)
!
*  ..locals
      logical :: lord
      logical, save :: lfirst=.true.
      integer :: i, j, ix, iy, ner, jak, jtemp, ipt, nxp2, n
      integer :: ixl, iyl, ixr, iyr, ixb, iyb, ixt, iyt
      integer :: ip(1:(nx+2)*(ny+2)+1)
      real (kind=R8) :: temp
      integer, save :: ncall = 0
      integer, save :: kwrite_fmt = -1
      integer, save :: kwrite_bnr = -1
      !$OMP THREADPRIVATE(lfirst, ncall, kwrite_fmt, kwrite_bnr)

      if(ncall.eq.0) then
        call ipgeti ('storage_11_kwrite_fmt', kwrite_fmt)
        call ipgeti ('storage_11_kwrite_bnr', kwrite_bnr)
        if (in_parallel() .and.
     &     (kwrite_bnr.ne.-1 .or. kwrite_fmt.ne.-1)) then
         !$OMP MASTER
         write(*,*) 'storage_11 OpenMP warning: no file IO ' ,
     &       'in parallel mode'
         !$OMP END MASTER
         kwrite_fmt = -1
         kwrite_bnr = -1
       endif
      endif

      n = (nx+2)*(ny+2)

      if (lfirst) then
        indexes=0
        nxp2=nx+2
        k=0
        ip(1)=1
        i=0
        do iy=-1,ny
          do ix=-1,nx
            ner=0
            i=i+1
*   ..Case of left neighbor
            jak = 2+leftix(ix,iy)+(2+leftiy(ix,iy)-1)*nxp2
            if( aa11(ix,iy,-1) .ne. 0.0e0_R8 ) then
              k = k+1
              a(k) = aa11(ix,iy,-1)
              indexes(1,k)=ix
              indexes(2,k)=iy
              indexes(3,k)=-1
              ia(k) = i
              ja(k) = jak
              ner = ner+1
            endif
*   ..Case (ix,iy)
            if( aa11(ix,iy,0) .ne. 0.0e0_R8 ) then
              k = k+1
              a(k) = aa11(ix,iy,0)
              indexes(1,k)=ix
              indexes(2,k)=iy
              indexes(3,k)=0
              ia(k) = i
              ja(k) = i
              ner = ner+1
            endif
*   ..Case of right neighbor
            jak = 2+rightix(ix,iy)+(2+rightiy(ix,iy)-1)*nxp2
            if( aa11(ix,iy,1) .ne. 0.0e0_R8 ) then
              k = k+1
              a(k) = aa11(ix,iy,1)
              indexes(1,k)=ix
              indexes(2,k)=iy
              indexes(3,k)=1
              ia(k) = i
              ja(k) = jak
              ner = ner+1
            endif
*   ..Case of top neighbor
            if( aa11(ix,iy,3) .ne. 0.0e0_R8 ) then
              k = k+1
              a(k) = aa11(ix,iy,3)
              indexes(1,k)=ix
              indexes(2,k)=iy
              indexes(3,k)=3
              ia(k) = i
              ja(k) = 2+topix(ix,iy)+(2+topiy(ix,iy)-1)*nxp2
              ner = ner+1
            endif
*   ..Case of bottom neighbor
            if( aa11(ix,iy,-3) .ne. 0.0e0_R8 ) then
              k = k+1
              a(k) = aa11(ix,iy,-3)
              indexes(1,k)=ix
              indexes(2,k)=iy
              indexes(3,k)=-3
              ia(k) = i
              ja(k) = 2+bottomix(ix,iy)+(2+bottomiy(ix,iy)-1)*nxp2
              ner = ner+1
            endif
*   ..Case (ix-1,iy-1)
            ixl = leftix(ix,iy)
            iyl = leftiy(ix,iy)
            ixb = bottomix(ix,iy)
            iyb = bottomiy(ix,iy)
            if( aa11(ix,iy,-4) .ne. 0.0e0_R8 ) then
              k = k+1
              a(k) = aa11(ix,iy,-4)
              indexes(1,k)=ix
              indexes(2,k)=iy
              indexes(3,k)=-4
              ia(k) = i
              ja(k) = 2+bottomix(ixl,iyl)+(2+bottomiy(ixl,iyl)-1)*nxp2
              ner = ner+1
            endif
            if( aa11(ix,iy,-6) .ne. 0.0e0_R8 ) then
              k = k+1
              a(k) = aa11(ix,iy,-6)
              indexes(1,k)=ix
              indexes(2,k)=iy
              indexes(3,k)= -6
              ia(k) = i
              ja(k) = 2+leftix(ixb,iyb)+(2+leftiy(ixb,iyb)-1)*nxp2
              ner = ner+1
            endif
*   ..Case (ix-1,iy+1)
            ixl = leftix(ix,iy)
            iyl = leftiy(ix,iy)
            ixt = topix(ix,iy)
            iyt = topiy(ix,iy)
            if( aa11(ix,iy,2) .ne. 0.0e0_R8 ) then
              k = k+1
              a(k) = aa11(ix,iy,2)
              indexes(1,k)=ix
              indexes(2,k)=iy
              indexes(3,k)=2
              ia(k) = i
              ja(k) = 2+topix(ixl,iyl)+(2+topiy(ixl,iyl)-1)*nxp2
              ner = ner+1
            endif
            if( aa11(ix,iy,5) .ne. 0.0e0_R8 ) then
              k = k+1
              a(k) = aa11(ix,iy,5)
              indexes(1,k)=ix
              indexes(2,k)=iy
              indexes(3,k)=5
              ia(k) = i
              ja(k) =2+leftix(ixt,iyt)+(2+leftiy(ixt,iyt)-1)*nxp2
              ner = ner+1
            endif
*       Case (ix+1,iy-1)
            ixr = rightix(ix,iy)
            iyr = rightiy(ix,iy)
            ixb = bottomix(ix,iy)
            iyb = bottomiy(ix,iy)
            if( aa11(ix,iy,-2) .ne. 0.0e0_R8 ) then
              k = k+1
              a(k) = aa11(ix,iy,-2)
              indexes(1,k)=ix
              indexes(2,k)=iy
              indexes(3,k)=-2
              ia(k) = i
              ja(k) =2+bottomix(ixr,iyr)+(2+bottomiy(ixr,iyr)-1)*nxp2
              ner = ner+1
            endif
            if( aa11(ix,iy,-5) .ne. 0.0e0_R8 ) then
              k = k+1
              a(k) = aa11(ix,iy,-5)
              indexes(1,k)=ix
              indexes(2,k)=iy
              indexes(3,k)=-5
              ia(k) = i
              ja(k) =2+rightix(ixb,iyb)+(2+rightiy(ixb,iyb)-1)*nxp2
              ner = ner+1
            endif
*       Case (ix+1,iy+1)
            ixr = rightix(ix,iy)
            iyr = rightiy(ix,iy)
            ixt = topix(ix,iy)
            iyt = topiy(ix,iy)
            if( aa11(ix,iy,4) .ne. 0.0e0_R8 ) then
              k = k+1
              a(k) = aa11(ix,iy,4)
              indexes(1,k)=ix
              indexes(2,k)=iy
              indexes(3,k)=4
              ia(k) = i
              ja(k) =2+topix(ixr,iyr)+(2+topiy(ixr,iyr)-1)*nxp2
              ner = ner+1
            endif
            if( aa11(ix,iy,6) .ne. 0.0e0_R8 ) then
              k = k+1
              a(k) = aa11(ix,iy,6)
              indexes(1,k)=ix
              indexes(2,k)=iy
              indexes(3,k)=6
              ia(k) = i
              ja(k) =2+rightix(ixt,iyt)+(2+rightiy(ixt,iyt)-1)*nxp2
              ner = ner+1
            endif
*   ..Case of top top neighbor
             ixt = topix(ix,iy)
             iyt = topiy(ix,iy)
            if( aa11(ix,iy,7) .ne. 0.0e0_R8 ) then
              k = k+1
              a(k) = aa11(ix,iy,7)
              indexes(1,k)=ix
              indexes(2,k)=iy
              indexes(3,k)=7
              ia(k) = i
              ja(k) =2+topix(ixt,iyt)+(2+topiy(ixt,iyt)-1)*nxp2
              ner = ner+1
            endif
*   ..Case of bottom bottom neighbor
             ixb = bottomix(ix,iy)
             iyb = bottomiy(ix,iy)
            if( aa11(ix,iy,-7) .ne. 0.0e0_R8 ) then
              k = k+1
              a(k) = aa11(ix,iy,-7)
              indexes(1,k)=ix
              indexes(2,k)=iy
              indexes(3,k)=-7
              ia(k) = i
              ja(k) = 2+bottomix(ixb,iyb)+(2+bottomiy(ixb,iyb)-1)*nxp2
              ner = ner+1
            endif
            ip(i+1)=ip(i) + ner

            do
              lord=.true.
              do ipt=ip(i),ip(i+1)-2
                if(ja(ipt).gt.ja(ipt+1)) then
                  jtemp  = ja(ipt)
                  ja(ipt)  = ja(ipt+1)
                  ja(ipt+1)= jtemp
                  temp   = a(ipt)
                  a(ipt)   = a(ipt+1)
                  a(ipt+1) = temp
                  do j=1,3
                    jtemp=indexes(j,ipt)
                    indexes(j,ipt)=indexes(j,ipt+1)
                    indexes(j,ipt+1)=jtemp
                  enddo
                  lord=.false.
                endif
              enddo
              if (lord) exit
            enddo

          enddo
        enddo

!       iasave=ia                                          !srv 15.05.07 {
!       jasave=ja
        do  i=1,k
          iasave(i)=ia(i)
          jasave(i)=ja(i)
        enddo                                              !srv 15.05.07 }
        lfirst=.false.

      else
!       ia=iasave                                           !srv 15.05.07 {
!       ja=jasave
        do i=1,k
          ia(i)=iasave(i)                                   !srv 15.05.07 }
          ja(i)=jasave(i)
          a(i)=aa11(indexes(1,i), indexes(2,i), indexes(3,i))
        enddo
      endif

      if(kwrite_bnr.eq.ncall) then
!        write(*,*) 'Opening 76'
!        open(76,file=
!     &   'st_11_rhs_binary_'//trim(name)//'.dat',form='unformatted')
!         write(76) n
!         do i=1,n
!          write(76) r1(i)
!         enddo
!        close(76)
        write(*,*) 'Opening 77'
        open(77,file=
     &   'st_11_matrix_binary_'//trim(name)//'.dat',form='unformatted') !srv 11.03.13
        write(77) n, k
        do i=1,k
          write(77) ia(i), ja(i), a(i)
        enddo
        close(77)
      endif
      if(kwrite_fmt.eq.ncall) then
!        write(*,*) 'Opening 76'
!        open(76,file='st_11_rhs_'//trim(name)//'.dat')
!         write(76,*) n
!         do i=1,n
!          write(76,*) r1(i)
!         enddo
!        close(76)
        open(77,file='st_11_matrix_'//trim(name)//'.dat')
        write(77,*) n, k
        do i=1,k
          write(77,*) ia(i), ja(i), a(i)
        enddo
        close(77)
      endif

      ncall = ncall + 1

      return
      end subroutine storage_11                                       !srv 02.04.18 }

