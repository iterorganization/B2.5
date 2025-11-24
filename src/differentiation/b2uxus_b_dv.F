*-----------------------------------------------------------------------
*.specification

      subroutine b2uxus_b_dv(nCv, mpg, aa, aad, aab, aabd, itcnt,   
     .  resfun, resfund, resfunb, resfunbd, corfun, corfund, corfunb,
     .  corfnubd, name, nbdirs)
      use b2mod_types
      use b2mod_ma28_for_us
      use b2us_map_diff_diffv
      use b2mod_lwmain
      use b2mod_subsys
      use b2mod_diffsizes
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, itcnt, nbdirs
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   aa(13*nCv), resfun(nCv)
      real (kind=R8) ::
     *   aab(13*nCv), resfunb(nCv)
      real (kind=R8) ::
     *   aad(nbdirsmax, 13*nCv), resfund(nbdirsmax, nCv)
      real (kind=R8) ::
     *   aabd(nbdirsmax, 13*nCv), resfunbd(nbdirsmax, nCv)
      character*(*) name 
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   corfun(nCv)
      real (kind=R8) ::
     *   corfunb(nCv)
      real (kind=R8) ::
     *   corfund(nbdirsmax, nCv)
      real (kind=R8) ::
     *   corfunbd(nbdirsmax, nCv)
*   ..common blocks

*-----------------------------------------------------------------------
*
*     Rearranged routine which solves a linear system of equation A*x=b
*     differentiated the first time in adjoint mode and the second in
*     tangent vector mode
*     The first adjoint derivation gives
*     AT = A'
*     AT*incrbb=xb
*     bb = bb + incrbb
*     A*x=b
*     Ab(i,j) = Ab(i,j) -x(j)*incrbb(i)
*     xb(i) = 0.0
*
*     The tangent differentiation gives
*     ATd = Ad'
*     AT = A'
*     d(AT*incrbb)=d(xb) --> SOLVE_DV(AT,ATd,xb,xbd,incrbb,incrbbd)
*     bbd = bbd + incrbbd
*     bb = bb + incrbb
*     d(A*x)=d(b) --> SOLVE_DV(A,Ad,b,bd,x,xd)
*     Abd(i,j) = Abd(i,j) -xd(j)*incrbb(i) -x(j)*incrbbd(i)
*     Ab(i,j) = Ab(i,j) -x(j)*incrbb(i)
*     xbd(i) = 0.0
*     xb(i) = 0.0
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, style, itmax, iters, itpar, cpu
      parameter (itmax=100)
      real (kind=R8) ::
     *   nrmc(0:itmax-1), cpu_start, cpu_end, c_end, c_start
      integer, save :: nz
      integer, allocatable, save :: ind(:), iasave(:), jasave(:)
      real (kind=R8) :: incrbb(nCv), incrbbd(nbdirsmax, nCv),
     &                  yy(nbdirsmax, nCv)
      integer, allocatable :: dummy(:)
      integer :: iCv, iNv
      integer mult_nonzero, mult_solvdim, mult_solvdim1
*   ..procedures
      intrinsic min, cpu_time
      external xertst, ipgeti, ipgetr,
     &  slv5pt, iluter, xerrab
*   ..initialization
      save ncall, style, cpu
      data ncall/0/, style/2/, cpu/0/
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2uxus_b_dv')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti ('b2ux9p_style', style)
       call xertst (0.le.style.and.style.lt.3,
     &   'faulty internal parameter style')
       call ipgeti ('b2ux9p_cpu', cpu)
       if(style.eq.2) then
         mult_nonzero=1
         mult_solvdim=15
         mult_solvdim1=10
         call ipgeti ('b2uxus_mult_nonzero', mult_nonzero)
         call ipgeti ('b2uxus_mult_solvdim', mult_solvdim)
         call ipgeti ('b2uxus_mult_solvdim1', mult_solvdim1)
         nnd=nCv
         n_nonzeros=mult_nonzero*mpg%nCmxNv
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
      call xertst (0.lt.nCv, 'faulty argument nCv')
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
     &                 ind, iasave, jasave, name)       
       call cpu_time(c_end)
*csc  Swap indexes to have transpose matrix
       dummy = irnh
       irnh = icnh
       icnh = dummy
       if(cpu.gt.0) write(*,'(a,f6.3)') 'cpu for storage_us ',
     &               c_end-c_start

*csc   Solve the system d(AT*bb)=d(xb) / d(AT*incrbb)=d(corfunb)
*      which gives ATd*bb + AT*bbd = xbd
*      so that AT*bbd = xbd - ATd*bb = yy
*      thus first solve AT*bb=xb  to obtain bb
       call ILUTERN_US(NCV, NZ, corfunb, incrbb, name)

*      Now do yy = xbd - ATd*bb for each direction
       do nd = 1, nbdirs
         yy(nd,:) = xbd(nd,:)
*        Call storage_us for matrix ATd
         call cpu_time(c_start)
         call storage_us (nCv, mpg, aad(nd,:), n_nonzeros, irnh, icnh,
     &                 a1h, nz, ind, iasave, jasave, name)       
         call cpu_time(c_end)
         if(cpu.gt.0) write(*,'(a,f6.3)') 'cpu for storage_us ',
     &               c_end-c_start
*        Swap indexes to have transpose matrix
         dummy = irnh
         irnh = icnh
         icnh = dummy
*        calculate yy = yy - ATd*bb
!      from Wouter
!      i,j,aij
!        j,b
!     do ii=1,nnz
         !get i,j
         ! ci = ci + aij*bj
         do nn = 1, n_nonzeros
           i = irnh(nn)
           j = icnh(nn)
           yy(nd,i) = yy(nd,i) - a1h(nn)*incrbb(j)
         end do
       end do
*      And finally solve AT*bbd = yy
*      Call storage_us for matrix ATd
       call cpu_time(c_start)
       call storage_us (nCv, mpg, aad(nd,:), n_nonzeros, irnh, icnh,
     &                 a1h, nz, ind, iasave, jasave, name)       
       call cpu_time(c_end)
       if(cpu.gt.0) write(*,'(a,f6.3)') 'cpu for storage_us ',
     &             c_end-c_start
*      Swap indexes to have transpose matrix
       dummy = irnh
       irnh = icnh
       icnh = dummy
       do nd = 1, nbdirs
        call ILUTERN_US(NCV, NZ, yy(nd,:), incrbbd(nd,:), name)
        resfunbd(nd,:) = resfunbd(nd,:) + incrbbd(nd,:)
       end do
       resfunb = resfunb + incrbb
      endif

*csc  solve d(A*x)=d(b) in tangent mode
      call b2uxus_dv (nCv, mpg, aa, aad, itcnt, resfun, resfund, corfun,
     & corfund, name, nbdirs)
*     Perform Abd(nd,i,j) = Abd(nd,i,j) -xd(nd,j)*incrbb(i)-x(j)*incrbbd(nd,i)
*     and Ab(i,j) = Ab(i,j) -x(j)*bb(i)
      do iCv=1,nCv
        do iNv=1,mpg%cvNvP(iCv,2)
          do nd=1, nbdirs
            aabd(nd,mpg%cvNvP(iCv,1)+iNv-1) = 
     &        aabd(nd,mpg%cvNvP(iCv,1)+iNv-1) - 
     &        corfund(nd,mpg%cvNv(mpg%cvNvP(iCv,1)+iNv-1))*incrbb(iCv) -
     &        corfun(mpg%cvNv(mpg%cvNvP(iCv,1)+iNv-1))*incrbbd(nd,iCv)
          end do
          aab(mpg%cvNvP(iCv,1)+iNv-1) = aab(mpg%cvNvP(iCv,1)+iNv-1) - 
     &      corfun(mpg%cvNv(mpg%cvNvP(iCv,1)+iNv-1))*incrbb(iCv)
        enddo
        corfunbd(1:nbdirs,iCv) = 0.0
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
*.end b2uxus_b_dv

      end subroutine b2uxus_b_dv

