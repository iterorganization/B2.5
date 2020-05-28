*-----------------------------------------------------------------------
*.specification

      subroutine b2ux7p (nx, ny, aa7, itcnt, resfun, corfun, name) !srv 11.01.00
      use b2mod_types
      use b2mod_ma28_for_7diag
      use b2mod_ysmp_sdrv
!      use b2mod_ma28
      use b2mod_indirect
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, itcnt
      real (kind=R8) ::
     *   aa7(-1:nx,-1:ny,-3:3), resfun(-1:nx,-1:ny)                !srv 11.01.00
      character*(*) name                                           !srv 18.05.02
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   corfun(-1:nx,-1:ny)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2UX7P solves a seven-point system of equations. This routine is
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
*     aa7 - (-1:nx,-1:ny,-3:3) real*8 array, input.           !srv 11.01.00
*     aa7 contains the seven-point matrix.                    !srv 11.01.00
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
      integer, save :: nz                                      !srv 29.05.02
      parameter (itmax=100)
      real (kind=R8) ::
     *   acpar, nrmc(0:itmax-1), cpu_start, cpu_end
      real (kind=R8)
     *  , allocatable, save :: aaf(:),wrk0(:),wrk1(:),wrk2(:)
      integer, allocatable, save :: iwrk(:)
      integer mult_nonzero, mult_solvdim, mult_solvdim1
*   ..procedures
      intrinsic min
      external subini, subend, xertst, ipgeti, ipgetr, slv5pt, iluter
*   ..initialisation
      save ncall, style, acpar, cpu
      data ncall/0/, style/2/, acpar/8.0e0_R8/, cpu/0/             !srv 26.12.17

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2ux7p')                                       !srv 11.01.00
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti ('b2ux7p_style', style)                         !srv 11.01.00
       call xertst (0.le.style.and.style.lt.4,
     &   'faulty internal parameter style')
       call ipgetr ('b2ux7p_acpar', acpar)                         !srv 11.01.00
       call xertst (0.lt.acpar, 'faulty internal parameter acpar')
       call ipgeti ('b2ux7p_cpu', cpu)                             !srv 11.01.00
       if(style.eq.2) then
         mult_nonzero=10
         mult_solvdim=15
         mult_solvdim1=10
         call ipgeti ('b2ux7p_mult_nonzero', mult_nonzero)         !srv 11.01.00
         call ipgeti ('b2ux7p_mult_solvdim', mult_solvdim)         !srv 11.01.00
         call ipgeti ('b2ux7p_mult_solvdim1', mult_solvdim1)       !srv 11.01.00
         nnd=(nx+2)*(ny+2)
         n_nonzeros=mult_nonzero*nnd
         solvdim=mult_solvdim*n_nonzeros
         solvdim1=mult_solvdim1*n_nonzeros
         call alloc_b2mod_ma28_for_7diag                           !srv 11.01.00
         call iniilu_7(nx+2,ny+2)                                  !srv 11.01.00
         nltrsol=2
         call ipgeti ('b2ux7p_nltrsol', nltrsol)                   !srv 11.01.00
       endif
       if(style.eq.3) then                                         !srv 29.05.02 {
         nnd=(nx+2)*(ny+2)
         nnd_sdrv  = (nx+2)*(ny+2)
         n_nz_sdrv = 4*nnd_sdrv
!         nsp = 3*nnd + 12*n_nz_sdrv  ! was 7*n_nz_sdrv without vis per !srv 13.02.03
         nsp = 3*nnd + 16*n_nz_sdrv                                !srv 29.12.15
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
     .   'Allocating space for solver for 7-point stencil matrix'  !srv 11.01.00
        if (style.eq.0) then
          call xerrab ('Style not applicable to 7-point stencil!')
*    ..obtain workspace for iluter
          allocate(aaf(1:n2*5))
          allocate(wrk0(1:n2))
          allocate(wrk1(1:n2))
          allocate(wrk2(1:n2))
          allocate(iwrk(1:1))
        else if (style.eq.1) then
          call xerrab ('Style not applicable to 7-point stencil!')
*    ..obtain workspace for slv5pt
          allocate(aaf(1:n2*(3*min(nx+2,ny+2)+1)))
          allocate(wrk0(1:n2))
          allocate(wrk1(1:1))
          allocate(wrk2(1:1))
          allocate(iwrk(1:n2))
        else if (style.eq.2) then
          continue
        else if (style.eq.3) then
          continue
        endif
        write(*,*) 
     .   'Allocated space for solver for 7-point stencil matrix'   !srv 11.01.00
      endif

*   ..solve the correction equation
      call cpu_time(cpu_start)
      if (style.eq.0) then
*    ..use iluter
       iters = min(itmax,4)
       itpar = itcnt*iters
       call iluter (nx+2, ny+2, aa7, resfun, iters, acpar, itpar,
     &   corfun, nrmc, aaf, wrk0, wrk1, wrk2)
      else if (style.eq.1) then
*    ..use slv5pt
       iters = 1
       call slv5pt (nx+2, ny+2, aa7, resfun, iters, corfun, nrmc,
     &   min(nx+2,ny+2), aaf, wrk0, iwrk)
      else if (style.eq.2) then
*    ..use ma28
        call ilutern_7(NX+2, NY+2, NX+2, aa7(-1,-1,-3),            !srv 11.01.00
     1   aa7(-1,-1,-2), aa7(-1,-1,-1), aa7(-1,-1,0),
     2   aa7(-1,-1,1), aa7(-1,-1,2), aa7(-1,-1,3),                 !srv 11.01.00
     .   resfun, corfun, name)
       else if (style.eq.3) then                                   !srv 29.05.02 {
!        call cpu_time(cpu_start)
        call storage_7 (nx, ny, aa7, n_nz_sdrv, ia, ja, a, nz,
     *                 indexes, iasave, jasave)
!        call cpu_time(cpu_end)
!        write (*,"(' Time of transform is',1pd9.2)") cpu_end-cpu_start
!        write(*,*) 'STORAGE_7: nz = ',nz
        call xertst (nz.le.N_NZ_SDRV, ' nz.gt.N_NZ_SDRV')
*    ..use sdrv from YSMP
!        call cpu_time(cpu_start)
        call illtern_7(NX+2, NY+2, NX+2, NZ, resfun, corfun)       !srv 29.05.02 }
!        call cpu_time(cpu_end)
!        write (*,"(' Time of illtern_7 is',1pd9.2)") cpu_end-cpu_start
      endif
      call cpu_time(cpu_end)
      if(cpu.gt.0) write(*,'(a,f6.3)') 'cpu for '//trim(name),
     &                                  cpu_end-cpu_start

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2ux7p

      end subroutine b2ux7p

      subroutine storage_7 (nx, ny, aa7, nnz, ip, ja, a, k,          !srv 29.05.02
     *                      indexes, ipsave, jasave)
      use b2mod_types
      use b2mod_indirect
      implicit none
*  ..input
      integer nx, ny, nnz
      real (kind=R8) :: aa7(-1:nx,-1:ny,-3:3)
*  ..output
      integer ip((nx+2)*(ny+2)+1), ja(nnz), k,
     *        indexes(1:3,nnz), ipsave((nx+2)*(ny+2)+1), jasave(nnz)
      real (kind=R8) :: a(nnz)
*  ..locals
      logical lord
      logical, save :: lfirst=.true.
      integer i, j, ix, iy, ner, iw, jw, jak, jtemp, ipt, nxp2
      real (kind=R8) :: temp

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
           jak = 2+leftix(ix,iy)+(2+leftiy(ix,iy)-1)*nxp2
           if( jak.gt.i .and. aa7(ix,iy,-1) .ne. 0.0e0_R8 ) then
             k = k+1
             a(k) = aa7(ix,iy,-1)
             indexes(1,k)=ix
             indexes(2,k)=iy
             indexes(3,k)=-1
             ja(k) = jak
             ner = ner+1
           endif
           if( aa7(ix,iy,0) .ne. 0.0e0_R8 ) then
             k = k+1
             a(k) = aa7(ix,iy,0)
             indexes(1,k)=ix
             indexes(2,k)=iy
             indexes(3,k)=0
             ja(k) = i
             ner = ner+1
           endif
           jak = 2+rightix(ix,iy)+(2+rightiy(ix,iy)-1)*nxp2
           if( jak.gt.i .and. aa7(ix,iy,1) .ne. 0.0e0_R8 ) then
             k = k+1
             a(k) = aa7(ix,iy,1)
             indexes(1,k)=ix
             indexes(2,k)=iy
             indexes(3,k)=1
             ja(k) = jak
             ner = ner+1
           endif
           if( aa7(ix,iy,2) .ne. 0.0e0_R8 ) then
             k = k+1
             a(k) = aa7(ix,iy,2)
             indexes(1,k)=ix
             indexes(2,k)=iy
             indexes(3,k)=2
             ja(k) =2+topix(ix,iy)+(2+topiy(ix,iy)-1)*nxp2
             ner = ner+1
           endif
           if( aa7(ix,iy,3) .ne. 0.0e0_R8 ) then
             k = k+1
             a(k) = aa7(ix,iy,3)
             indexes(1,k)=ix
             indexes(2,k)=iy
             indexes(3,k)=3
             iw = topix(ix,iy)
             jw = topiy(ix,iy)
             ja(k) =2+topix(iw,jw)+(2+topiy(iw,jw)-1)*nxp2
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

       ipsave=ip
       jasave=ja
       lfirst=.false.

      else
       ip=ipsave
       ja=jasave
       do i=1,k
        a(i)=aa7(indexes(1,i), indexes(2,i), indexes(3,i))
       enddo
      endif

      end                                                         !srv 29.05.02
