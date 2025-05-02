*-----------------------------------------------------------------------
*.specification

      subroutine b2usht (nCv, nFc, nVx, ns, switch, geo, mpg,
     &  itcnt, ne, ni, nn, te, ti, tn, kt, zt,
     &  floe, cone, floi, coni, flon, conn, flokt, conkt, flozt, conzt,
     &  she, shi, shn, skt, szt, ceqp,
     &  resht, reshe, reshi, reshn, reskt, reszt,
     &  cortt, corte, corti, cortn, corkt, corzt, flo0, con0, aa)
      use b2mod_types
      use b2mod_numerics_namelist
      use b2mod_switches
      use b2mod_constants
      use b2us_geo
      use b2us_map
      use b2mod_ma28_for_us
     , , only : restart_ma28_for_us
      use b2mod_openmp
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2usht
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, itcnt
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *  ne(nCv), ni(nCv,0:1), nn(nCv), te(nCv), ti(nCv), tn(nCv),
     &  kt(nCv), zt(nCv), floe(nFc,0:1), floi(nFc,0:1), flon(nFc,0:1),
     &  flokt(nFc,0:1), flozt(nFc,0:1),
     &  cone(nFc,0:1,0:2), coni(nFc,0:1,0:2),
     &  conn(nFc,0:1), conkt(nFc,0:1), conzt(nFc,0:1)
      real (kind=R8) ::
     *  she(nCv,0:3), shi(nCv,0:3), shn(nCv,0:3), skt(nCv,0:3),
     &  szt(nCv,0:3),
     &  ceqp(nCv), resht(nCv), reshe(nCv), reshi(nCv), reshn(nCv),
     &  reskt(nCv), reszt(nCv)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  cortt(nCv), corte(nCv), corti(nCv), cortn(nCv), corkt(nCv),
     &  corzt(nCv)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *  flo0(nFc,0:1), con0(nFc,0:1,0:2), aa(mpg%nCmxNv)
      real (kind=R8) :: dumm(1)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2USHT computes a correction to the fields te and ti.
*
*-----------------------------------------------------------------------
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
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, iFc, iCv1, iCv2, iVx1, iVx2, inCv
c      integer, save :: no_Ptncr_x_he = 0, no_Ptncr_y_he = 0            !srv 06.07.06
c      integer, save :: no_Ptncr_x_hi = 0, no_Ptncr_y_hi = 0            !srv 06.07.06
c      integer, save :: iadd_9 = 1                                      !srv 02.01.07
      integer ireg
      logical lnonzero
      real (kind=R8) ::
     *   t0, t1, t2, t3, t11, t22
      logical new_matrix
#ifdef _OPENMP
      logical make_parallel
#endif 
*   ..procedures
      intrinsic abs, sqrt, dim
      external xertst, sfill
      external b2xvsg, b2xvff, b2xvfx, b2uxus                           !srv 02.01.07
*   ..initialization

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
#ifdef _OPENMP
      ! Check if we're already running in parallel 
      make_parallel = in_parallel()
#endif 
      if (switch%b2usht_iout.ne.0) then  
        if (in_parallel() ) then 
          if (ncall_b2usht == 0) then 
#ifdef _OPENMP
!$OMP MASTER
#endif
      write(*,*) 'b2ush9 OpenMP warning: no file IO in parallel mode'
#ifdef _OPENMP
!$OMP END MASTER
#endif
          endif 
        endif 
      endif
      iCv1 = 1 !csc needed for correct initialization in adjoint AD
      iCv2 = 1
*   ..subprogram start-up calls
      call subini ('b2usht')
*   ..test nCv, nFc, ns
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test itcnt, rxg
      call xertst (0.le.itcnt, 'faulty argument itcnt')
      call xertst (0.0_R8.lt.switch%b2npht_rxg, 'faulty argument rxg')
*   ..extensive tests on first few calls
      if (ncall_b2usht.lt.3) then
*    ..test sign of ne, ni, te, ti
       call b2xvsg (nCv, ne, 1, 'ne', '.gt.')
       call b2xvsg (nCv, ni, 1, 'ni', '.gt.')
       call b2xvsg (nCv, te, 1, 'te', '.gt.')
       call b2xvsg (nCv, ti, 1, 'ti', '.gt.')
       call b2xvsg (nCv, kt, 1, 'kt', '.ge.')
       call b2xvsg (nCv, zt, 1, 'zt', '.ge.')
*    ..test sign of she, shi, ceqp
       if (switch%b2mndt_style.ne.2) then
         call b2xvsg (nCv, she(1,0), 1, 'she0', '.ge.')
       endif
       call b2xvsg (nCv, she(1,1), 1, 'she1', '.le.')
       call b2xvsg (nCv, she(1,2), 1, 'she2', '.ge.')
       call b2xvsg (nCv, she(1,3), 1, 'she3', '.le.')
       if (switch%b2mndt_style.ne.2) then
         call b2xvsg (nCv, shi(1,0), 1, 'shi0', '.ge.')
       endif
       call b2xvsg (nCv, shi(1,1), 1, 'shi1', '.le.')
       call b2xvsg (nCv, shi(1,2), 1, 'shi2', '.ge.')
       call b2xvsg (nCv, shi(1,3), 1, 'shi3', '.le.')
       if (switch%b2mndt_style.ne.2) then
         call b2xvsg (nCv, shn(1,0), 1, 'shn0', '.ge.')
       endif
       call b2xvsg (nCv, shn(1,1), 1, 'shn1', '.le.')
       call b2xvsg (nCv, shn(1,2), 1, 'shn2', '.ge.')
       call b2xvsg (nCv, shn(1,3), 1, 'shn3', '.le.')
       call b2xvsg (nCv, ceqp, 1, 'ceqp', '.ge.')
      endif

      if (switch%tn_style.le.1) then !! in case we solve a combined ion-neutral energy equation (default) or pure ion energy equation we make sure
         conn = 0.0_R8               !! the neutral fluxes and sources are zero.
         flon = 0.0_R8
         shn = 0.0_R8
      endif

#if (defined(_OPENMP) && !defined(NO_OPENMP_B2USHT))
* ..initialize openmp parallel region for workshare
!svdk the current implementation permits to parallellize this section
!svdk using a task construct, since each energy equation is solved
!svdk separately first, and updates are applied only later.
!$omp parallel default(none) 
!$omp& shared(mpg, switch, geo, nCv, nFc, itcnt)
!$omp& shared(solveee, solveei, solveen, solveet, solvekt, solvezt)
!$omp& shared(floe, floi, flon, flokt, flozt)
!$omp& shared(cone, coni, conn, conkt, conzt)
!$omp& shared(resht, cortt, reshe, corte, reshi, corti, reshn, cortn)
!$omp& shared(reskt, corkt, reszt, corzt)
!$omp& shared(she, shi, shn, ne, ni, nn, te, ti, tn, ceqp)
!$omp& shared(skt, kt, szt, zt)
!$omp& shared(last_solve_9) 
!$omp& private(aa, flo0, con0)
!$omp& private(iCv, iFc, iCv1, iCv2, iVx1, iVx2, inCv)
!$omp& private(ireg, t0, t1, t2, t11, t22, t3, lnonzero, dumm)
!svdk task 1: total heat correction equation
!svdk should be ran by a single thread, tasks are accumulated and 
!svdk executed by other threads
!!$omp single 
!!$omp task
!$omp sections
!$omp section
#endif 
* ..total heat correction equation 
      if (any(solveee(0:mpg%nnreg(0))).and.
     &    any(solveei(0:mpg%nnreg(0))).and.
     &    any(solveen(0:mpg%nnreg(0))).and.
     &    any(solveet(0:mpg%nnreg(0)))) then                                !srv 22.05.18
   
#ifdef _OPENMP 
      call restart_ma28_for_us
#else
      new_matrix = .false.
      do ireg = 0, mpg%nnreg(0)
        new_matrix = new_matrix.or.
     &   solveet(ireg).neqv.last_solve_9(ireg)
      end do
      if (new_matrix.or..true.) call restart_ma28_for_us
#endif
* ..total heat correction equation
*   ..compute flo0, con0
      flo0        = floe + floi + flon
      con0(:,:,0) = cone(:,:,0) + coni(:,:,0)+ conn
      if (switch%b2usht_iout.ne.0) then  
          if (ncall_b2usht == 0 .and. .not. in_parallel()) then 
            call my_out_us(70,nFc,1,flo0(1,0),'b2usht_flo0x')
            call my_out_us(70,nFc,1,flo0(1,1),'b2usht_flo0y')
            call my_out_us(70,nFc,1,con0(1,0,0),'b2usht_con0xtest')
            call my_out_us(70,nFc,1,con0(1,1,0),'b2usht_con0ytest')
            call my_out_us(70,nFc,1,con0(1,0,1),'b2usht_con0x10')
            call my_out_us(70,nFc,1,con0(1,1,1),'b2usht_con0y10')
            call my_out_us(70,nFc,1,con0(1,0,2),'b2usht_con0x20') 
            call my_out_us(70,nFc,1,con0(1,1,2),'b2usht_con0y20') 
        end if 
      endif

*   ..initialize matrix to zero
      aa = 0.0_R8
*   ..compute the nine-point matrix
! $AD II-LOOP
      do iCv = 1, nCv
        ireg = mpg%cvReg(iCv)
        t0 = 1.0e60_R8
        if (solveee(ireg)) then
          aa(mpg%cvNvP(iCv,1)) = aa(mpg%cvNvP(iCv,1))
     &    -she(iCv,1)-she(iCv,3)*ne(iCv)
          t0 = min(t0,te(iCv))
        endif
        if (solveei(ireg)) then
          aa(mpg%cvNvP(iCv,1)) = aa(mpg%cvNvP(iCv,1))
     &    -shi(iCv,1)-shi(iCv,3)*ni(iCv,0)
          t0 = min(t0,ti(iCv))
        endif
        if (solveen(ireg).and.switch%tn_style.eq.2) then
          aa(mpg%cvNvP(iCv,1)) = aa(mpg%cvNvP(iCv,1))
     &    -shn(iCv,1)-shn(iCv,3)*nn(iCv)
          t0 = min(t0,tn(iCv))
        endif
        if (solveee(ireg).or.solveei(ireg).or.
     &      (solveen(ireg).and.switch%tn_style.eq.2)) then
          if (switch%b2mndt_style.lt.2) then
            aa(mpg%cvNvP(iCv,1)) = aa(mpg%cvNvP(iCv,1))
     &        +abs(resht(iCv))/(switch%b2npht_rxg*t0)
          endif
        endif
      enddo
! $AD II-LOOP

      if (switch%b2usht_iout.ne.0 .and. .not. in_parallel()) then
       call my_out_us(70,nCv,0,she(1,1),'b2usht_she1')     
       call my_out_us(70,nCv,0,she(1,3),'b2usht_she3')
       call my_out_us(70,nCv,0,shi(1,1),'b2usht_shi1')     
       call my_out_us(70,nCv,0,shi(1,3),'b2usht_shi3')
       call my_out_us(70,nCv,0,shn(1,1),'b2usht_shn1')     
       call my_out_us(70,nCv,0,shn(1,3),'b2usht_shn3')
       call my_out_us(70,nCv,0,ne,'b2usht_ne')
       call my_out_us(70,nCv,0,ni,'b2usht_ni')
       call my_out_us(70,nCv,0,te,'b2usht_te')
       call my_out_us(70,nCv,0,ti,'b2usht_ti')
       call my_out_us(70,nCv,0,resht,'b2usht_resht')

       call my_out_us(70,mpg%nCmxNv,0,aa,'b2usht_aa_tt_pre')
      endif

      do iFc = 1, nFc
        iCv1 = mpg%fcCv(iFc,1) ! cell number of c1
        iCv2 = mpg%fcCv(iFc,2) ! cell number of c2
*     ..contributions from flo0(:,0), con0(:,0,:) ('5-point' contributions)
        t1 = - 0.5_R8*flo0(iFc,0) - con0(iFc,0,0)
        t2 =   0.5_R8*flo0(iFc,0) - con0(iFc,0,0)
        if (solveee(mpg%cvReg(iCv1)).or.solveei(mpg%cvReg(iCv1)).or.
     &      (solveen(ireg).and.switch%tn_style.eq.2)) then
         aa(mpg%cvNvP(iCv1,1)) = aa(mpg%cvNvP(iCv1,1)) - t1
         where (mpg%cvNv(mpg%cvNvP(iCv1,1):
     &                   mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) == iCv2)
          aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) =
     &     aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) +
     &     t2
         end where
        endif
        if (solveee(mpg%cvReg(iCv2)).or.solveei(mpg%cvReg(iCv2)).or.
     &      (solveen(ireg).and.switch%tn_style.eq.2)) then
         aa(mpg%cvNvP(iCv2,1)) = aa(mpg%cvNvP(iCv2,1)) - t2
         where (mpg%cvNv(mpg%cvNvP(iCv2,1):
     &                   mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) == iCv1)
          aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) =
     &     aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) +
     &     t1
         end where
        endif
*     ..contributions from flo0(:,1), con0(:,1,:) ('9-point' contributions)
        iVx1 = mpg%fcVx(iFc,1) ! number of v1
        iVx2 = mpg%fcVx(iFc,2) ! number of v2
        t1 = - 0.5_R8*flo0(iFc,1) - con0(iFc,1,0)
        t2 =   0.5_R8*flo0(iFc,1) - con0(iFc,1,0)
        t11 = 0.0_R8
        do inCv = 1, mpg%vxCvP(iVx1,2)
          t11 = t11 + 1/geo%vxVol(mpg%vxCvP(iVx1,1)+inCv-1)
        enddo
        t22 = 0.0_R8
        do inCv = 1, mpg%vxCvP(iVx2,2)
          t22 = t22 + 1/geo%vxVol(mpg%vxCvP(iVx2,1)+inCv-1)
        enddo
        ! contribution to equation of iCv1
        if (solveee(mpg%cvReg(iCv1)).or.solveei(mpg%cvReg(iCv1)).or.
     &      (solveen(ireg).and.switch%tn_style.eq.2)) then
         ! contribution from iVx1
         do inCv = 1, mpg%vxCvP(iVx1,2)
          iCv = mpg%vxCv(mpg%vxCvP(iVx1,1) + inCv - 1)
          t3  = 1/(geo%vxVol(mpg%vxCvP(iVx1,1) + inCv - 1)*t11)
          where (mpg%cvNv(mpg%cvNvP(iCv1,1):
     &                    mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) == iCv)
           aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) =
     &      aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1)-
     &      t3*t1
          end where
         enddo
         ! contribution from iVx2
         do inCv = 1, mpg%vxCvP(iVx2,2)
          iCv = mpg%vxCv(mpg%vxCvP(iVx2,1) + inCv - 1)
          t3  = 1/(geo%vxVol(mpg%vxCvP(iVx2,1) + inCv - 1)*t22)
          where (mpg%cvNv(mpg%cvNvP(iCv1,1):
     &                    mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) == iCv)
           aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) =
     &      aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1)+
     &      t3*t2
          end where
         enddo
        endif
        ! contribution to equation of iCv2
        if (solveee(mpg%cvReg(iCv2)).or.solveei(mpg%cvReg(iCv2)).or.
     &      (solveen(ireg).and.switch%tn_style.eq.2)) then
         ! contribution from iVx1
         do inCv = 1, mpg%vxCvP(iVx1,2)
          iCv = mpg%vxCv(mpg%vxCvP(iVx1,1) + inCv - 1)
          t3  = 1/(geo%vxVol(mpg%vxCvP(iVx1,1) + inCv - 1)*t11)
          where (mpg%cvNv(mpg%cvNvP(iCv2,1):
     &                    mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) == iCv)
           aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) =
     &      aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1)+
     &      t3*t1
          end where
         enddo
         ! contribution from iVx2
         do inCv = 1, mpg%vxCvP(iVx2,2)
          iCv = mpg%vxCv(mpg%vxCvP(iVx2,1) + inCv - 1)
          t3  = 1/(geo%vxVol(mpg%vxCvP(iVx2,1) + inCv - 1)*t22)
          where (mpg%cvNv(mpg%cvNvP(iCv2,1):
     &                    mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) == iCv)
           aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) =
     &      aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1)-
     &      t3*t2
          end where
         enddo
        endif
      enddo
*   ..check for empty matrix lines and replace with dummy equation
! $AD II-LOOP
      do iCv = 1, nCv
        lnonzero=.false.
        do inCv = 1,mpg%cvNvP(iCv,2)
          lnonzero = lnonzero.or.(aa(mpg%cvNvP(iCv,1)+inCv-1).ne.0.0_R8)
        enddo
        if (.not.lnonzero) then
          aa(mpg%cvNvP(iCv,1)) = 1.0_R8
          resht(iCv) = 0.0_R8
        endif
      enddo
      if (switch%b2usht_iout.ne.0 .and. .not. in_parallel()) then
        call my_out_us(70,mpg%nCmxNv,0,aa,'b2usht_aa_tt')
      endif 
*   ..solve the correction equation
      call b2uxus (nCv, mpg, aa, itcnt, resht, cortt, 'b2usht_total')  !srv 18.05.02
      endif                                                             !srv 22.05.18

      if (switch%b2usht_iout.ne.0 .and. .not. in_parallel()) then
       call my_out_us(70,mpg%nCmxNv,0,aa,'b2usht_ht_aa')
       call my_out_us(70,nCv,0,resht,'b2usht_resht')
       call my_out_us(70,nCv,0,cortt,'b2usht_cortt')
      endif

#if (defined(_OPENMP) && !defined(NO_OPENMP_B2USHT))
!!$omp end task
!svdk task 2: electron temperature equation
!!$omp task
!$omp section
#endif 
* ..electron temperature correction equation
      if (any(solveee(0:mpg%nnreg(0)))) then                                 !srv 22.05.18
#ifdef _OPENMP 
      call restart_ma28_for_us
#else 
      new_matrix = .false.
      do ireg = 0, mpg%nnreg(0)
        new_matrix = new_matrix.or.
     &   solveee(ireg).neqv.solveet(ireg)
      end do
      if (new_matrix.or..true.) call restart_ma28_for_us
#endif 
*   ..compute flo0, con0
      flo0        = floe
      con0(:,:,0) = cone(:,:,0)
*   ..initialize matrix to zero
      aa = 0.0_R8
*   ..compute the nine-point matrix
! $AD II-LOOP
      do iCv = 1, nCv
        ireg = mpg%cvReg(iCv)
        if (solveee(ireg)) then
          aa(mpg%cvNvP(iCv,1)) =
     &    -she(iCv,1)-she(iCv,3)*ne(iCv)+
     &    dim(1.0_R8,1.5_R8*(te(iCv)-ti(iCv))/te(iCv))*               !srv 02.01.07
     &    ceqp(iCv)*ne(iCv)*ni(iCv,0)/(te(iCv)*sqrt(te(iCv)))
          if (switch%b2mndt_style.lt.2) then
            aa(mpg%cvNvP(iCv,1)) = aa(mpg%cvNvP(iCv,1))
     &      +abs(reshe(iCv))/(switch%b2npht_rxg*te(iCv))              !srv 02.01.07
          endif
        endif
      enddo
! $AD II-LOOP

      if (switch%b2usht_iout.ne.0 .and. .not. in_parallel()) then
       call my_out_us(70,nCv,0,she(1,1),'b2usht_she1_te')     
       call my_out_us(70,nCv,0,she(1,3),'b2usht_she3_te')
       call my_out_us(70,nCv,0,te,'b2usht_te_te')
       call my_out_us(70,nCv,0,ti,'b2usht_ti_te')
       call my_out_us(70,nCv,0,ne,'b2usht_ne_te')
       call my_out_us(70,nCv,0,ni,'b2usht_ni_te')
       call my_out_us(70,nCv,0,reshe,'b2usht_reshe_te')
       dumm = real(mpg%nCmxNv,R8)
       call my_out_us(70,1,0,dumm,'b2usht_nCmxNv')

       call my_out_us(70,nFc,1,flo0(1,0),'b2usht_floex') 
       call my_out_us(70,nFc,1,flo0(1,1),'b2usht_floey')  
       call my_out_us(70,nFc,1,con0(1,0,0),'b2usht_conex')
       call my_out_us(70,nFc,1,con0(1,1,0),'b2usht_coney')
       call my_out_us(70,mpg%nCmxNv,0,aa,'b2usht_aa_te_pre')
      endif 

      do iFc = 1, nFc
        iCv1 = mpg%fcCv(iFc,1) ! cell number of c1
        iCv2 = mpg%fcCv(iFc,2) ! cell number of c2
*     ..contributions from flo0(:,0), con0(:,0,:) ('5-point' contributions)
        t1 = - 0.5_R8*flo0(iFc,0) - con0(iFc,0,0)
        t2 =   0.5_R8*flo0(iFc,0) - con0(iFc,0,0)
        if (solveee(mpg%cvReg(iCv1))) then
         aa(mpg%cvNvP(iCv1,1)) = aa(mpg%cvNvP(iCv1,1)) - t1
         where (mpg%cvNv(mpg%cvNvP(iCv1,1):
     &                   mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) == iCv2)
          aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) =
     &     aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) +
     &     t2
         end where
        endif
        if (solveee(mpg%cvReg(iCv2))) then
         aa(mpg%cvNvP(iCv2,1)) = aa(mpg%cvNvP(iCv2,1)) - t2
         where (mpg%cvNv(mpg%cvNvP(iCv2,1):
     &                   mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) == iCv1)
          aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) =
     &     aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) +
     &     t1
         end where
        endif
*     ..contributions from flo0(:,1), con0(:,1,:) ('9-point' contributions)
        iVx1 = mpg%fcVx(iFc,1) ! number of v1
        iVx2 = mpg%fcVx(iFc,2) ! number of v2
        t1 = - 0.5_R8*flo0(iFc,1) - con0(iFc,1,0)
        t2 =   0.5_R8*flo0(iFc,1) - con0(iFc,1,0)
        t11 = 0.0_R8
        do inCv = 1, mpg%vxCvP(iVx1,2)
          t11 = t11 + 1/geo%vxVol(mpg%vxCvP(iVx1,1)+inCv-1)
        enddo
        t22 = 0.0_R8
        do inCv = 1, mpg%vxCvP(iVx2,2)
          t22 = t22 + 1/geo%vxVol(mpg%vxCvP(iVx2,1)+inCv-1)
        enddo
        ! contribution to equation of iCv1
        if (solveee(mpg%cvReg(iCv1))) then
         ! contribution from iVx1
         do inCv = 1, mpg%vxCvP(iVx1,2)
          iCv = mpg%vxCv(mpg%vxCvP(iVx1,1) + inCv - 1)
          t3  = 1/(geo%vxVol(mpg%vxCvP(iVx1,1) + inCv - 1)*t11)
          where (mpg%cvNv(mpg%cvNvP(iCv1,1):
     &                    mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) == iCv)
           aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) =
     &      aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1)-
     &      t3*t1
          end where
         enddo
         ! contribution from iVx2
         do inCv = 1, mpg%vxCvP(iVx2,2)
          iCv = mpg%vxCv(mpg%vxCvP(iVx2,1) + inCv - 1)
          t3  = 1/(geo%vxVol(mpg%vxCvP(iVx2,1) + inCv - 1)*t22)
          where (mpg%cvNv(mpg%cvNvP(iCv1,1):
     &                    mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) == iCv)
           aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) =
     &      aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1)+
     &      t3*t2
          end where
         enddo
        endif
        ! contribution to equation of iCv2
        if (solveee(mpg%cvReg(iCv2))) then
         ! contribution from iVx1
         do inCv = 1, mpg%vxCvP(iVx1,2)
          iCv = mpg%vxCv(mpg%vxCvP(iVx1,1) + inCv - 1)
          t3  = 1/(geo%vxVol(mpg%vxCvP(iVx1,1) + inCv - 1)*t11)
          where (mpg%cvNv(mpg%cvNvP(iCv2,1):
     &                    mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) == iCv)
           aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) =
     &      aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1)+
     &      t3*t1
          end where
         enddo
         ! contribution from iVx2
         do inCv = 1, mpg%vxCvP(iVx2,2)
          iCv = mpg%vxCv(mpg%vxCvP(iVx2,1) + inCv - 1)
          t3  = 1/(geo%vxVol(mpg%vxCvP(iVx2,1) + inCv - 1)*t22)
          where (mpg%cvNv(mpg%cvNvP(iCv2,1):
     &                    mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) == iCv)
           aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) =
     &      aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1)-
     &      t3*t2
          end where
         enddo
        endif
      enddo
*   ..check for empty matrix lines and replace with dummy equation
! $AD II-LOOP
      do iCv = 1, nCv
        lnonzero=.false.
        do inCv = 1,mpg%cvNvP(iCv,2)
          lnonzero = lnonzero.or.(aa(mpg%cvNvP(iCv,1)+inCv-1).ne.0.0_R8)
        enddo
        if (.not.lnonzero) then
          aa(mpg%cvNvP(iCv,1)) = 1.0_R8
          reshe(iCv) = 0.0_R8
        endif
      enddo

      if (switch%b2usht_iout.ne.0 .and. .not. in_parallel()) then
       call my_out_us(70,mpg%nCmxNv,0,aa,'b2usht_aa_te')
      endif
*   ..solve the correction equation
      call b2uxus (nCv, mpg, aa, itcnt, reshe, corte, 'b2usht_Te')     !srv 18.05.02
      endif                                                            !srv 22.05.18

      if (switch%b2usht_iout.ne.0 .and. .not. in_parallel()) then
       call my_out_us(70,mpg%nCmxNv,0,aa,'b2usht_he_aa')
       call my_out_us(70,nCv,0,reshe,'b2usht_reshe')
       call my_out_us(70,nCv,0,corte,'b2usht_corte')
      endif

#if (defined(_OPENMP) && !defined(NO_OPENMP_B2USHT))
!!$omp end task
!svdk task 3: ion temperature equation
!!$omp task
!$omp section
#endif 
* ..atom temperature correction equation
      if (any(solveei(0:mpg%nnreg(0)))) then                                 !srv 22.05.18
#ifdef _OPENMP 
      call restart_ma28_for_us  
#else 
      new_matrix = .false.
      do ireg = 0, mpg%nnreg(0)
        new_matrix = new_matrix.or.
     &   solveei(ireg).neqv.solveee(ireg)
      end do
      if (new_matrix.or..true.) call restart_ma28_for_us
#endif 
*   ..compute flo0, con0
      flo0        = floi
      con0(:,:,0) = coni(:,:,0)
*   ..initialize matrix to zero
      aa = 0.0_R8
*   ..compute the nine-point matrix
! $AD II-LOOP
      do iCv = 1, nCv
        ireg = mpg%cvReg(iCv)
        if (solveei(ireg)) then
          aa(mpg%cvNvP(iCv,1)) =
     &    -shi(iCv,1)-shi(iCv,3)*ni(iCv,0)+
     &    ceqp(iCv)*ne(iCv)*ni(iCv,0)/(te(iCv)*sqrt(te(iCv)))
          if (switch%b2mndt_style.lt.2) then
            aa(mpg%cvNvP(iCv,1)) = aa(mpg%cvNvP(iCv,1))
     &      +abs(reshi(iCv))/(switch%b2npht_rxg*ti(iCv))
          endif
        endif
      enddo
! $AD II-LOOP

      if (switch%b2usht_iout.ne.0 .and. .not. in_parallel()) then
       call my_out_us(70,nFc,1,floi(1,0),'b2usht_floix') 
       call my_out_us(70,nFc,1,floi(1,1),'b2usht_floiy')  
       call my_out_us(70,nFc,1,coni(1,0,0),'b2usht_conix')
       call my_out_us(70,nFc,1,coni(1,1,0),'b2usht_coniy')

       call my_out_us(70,nCv,0,shi(1,1),'b2usht_shi1_ti')     
       call my_out_us(70,nCv,0,shi(1,3),'b2usht_shi3_ti')
       call my_out_us(70,nCv,0,te,'b2usht_te_ti')
       call my_out_us(70,nCv,0,ti,'b2usht_ti_ti')
       call my_out_us(70,nCv,0,ne,'b2usht_ne_ti')
       call my_out_us(70,nCv,0,ni,'b2usht_ni_ti')
       call my_out_us(70,nCv,0,reshi,'b2usht_reshi_ti')

       call my_out_us(70,mpg%nCmxNv,0,aa,'b2usht_aa_ti_pre')
      endif

      do iFc = 1, nFc
        iCv1 = mpg%fcCv(iFc,1) ! cell number of c1
        iCv2 = mpg%fcCv(iFc,2) ! cell number of c2
*     ..contributions from flo0(:,0), con0(:,0,:) ('5-point' contributions)
        t1 = - 0.5_R8*flo0(iFc,0) - con0(iFc,0,0)
        t2 =   0.5_R8*flo0(iFc,0) - con0(iFc,0,0)
        if (solveei(mpg%cvReg(iCv1))) then
         aa(mpg%cvNvP(iCv1,1)) = aa(mpg%cvNvP(iCv1,1)) - t1
         where (mpg%cvNv(mpg%cvNvP(iCv1,1):
     &                   mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) == iCv2)
          aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) =
     &     aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) +
     &     t2
         end where
        endif
        if (solveei(mpg%cvReg(iCv2))) then
         aa(mpg%cvNvP(iCv2,1)) = aa(mpg%cvNvP(iCv2,1)) - t2
         where (mpg%cvNv(mpg%cvNvP(iCv2,1):
     &                   mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) == iCv1)
          aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) =
     &     aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) +
     &     t1
         end where
        endif
*     ..contributions from flo0(:,1), con0(:,1,:) ('9-point' contributions)
        iVx1 = mpg%fcVx(iFc,1) ! number of v1
        iVx2 = mpg%fcVx(iFc,2) ! number of v2
        t1 = - 0.5_R8*flo0(iFc,1) - con0(iFc,1,0)
        t2 =   0.5_R8*flo0(iFc,1) - con0(iFc,1,0)
        t11 = 0.0_R8
        do inCv = 1, mpg%vxCvP(iVx1,2)
          t11 = t11 + 1/geo%vxVol(mpg%vxCvP(iVx1,1)+inCv-1)
        enddo
        t22 = 0.0_R8
        do inCv = 1, mpg%vxCvP(iVx2,2)
          t22 = t22 + 1/geo%vxVol(mpg%vxCvP(iVx2,1)+inCv-1)
        enddo
        ! contribution to equation of iCv1
        if (solveei(mpg%cvReg(iCv1))) then
         ! contribution from iVx1
         do inCv = 1, mpg%vxCvP(iVx1,2)
          iCv = mpg%vxCv(mpg%vxCvP(iVx1,1) + inCv - 1)
          t3  = 1/(geo%vxVol(mpg%vxCvP(iVx1,1) + inCv - 1)*t11)
          where (mpg%cvNv(mpg%cvNvP(iCv1,1):
     &                    mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) == iCv)
           aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) =
     &      aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1)-
     &      t3*t1
          end where
         enddo
         ! contribution from iVx2
         do inCv = 1, mpg%vxCvP(iVx2,2)
          iCv = mpg%vxCv(mpg%vxCvP(iVx2,1) + inCv - 1)
          t3  = 1/(geo%vxVol(mpg%vxCvP(iVx2,1) + inCv - 1)*t22)
          where (mpg%cvNv(mpg%cvNvP(iCv1,1):
     &                    mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) == iCv)
           aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) =
     &      aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1)+
     &      t3*t2
          end where
         enddo
        endif
        ! contribution to equation of iCv2
        if (solveei(mpg%cvReg(iCv2))) then
         ! contribution from iVx1
         do inCv = 1, mpg%vxCvP(iVx1,2)
          iCv = mpg%vxCv(mpg%vxCvP(iVx1,1) + inCv - 1)
          t3  = 1/(geo%vxVol(mpg%vxCvP(iVx1,1) + inCv - 1)*t11)
          where (mpg%cvNv(mpg%cvNvP(iCv2,1):
     &                    mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) == iCv)
           aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) =
     &      aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1)+
     &      t3*t1
          end where
         enddo
         ! contribution from iVx2
         do inCv = 1, mpg%vxCvP(iVx2,2)
          iCv = mpg%vxCv(mpg%vxCvP(iVx2,1) + inCv - 1)
          t3  = 1/(geo%vxVol(mpg%vxCvP(iVx2,1) + inCv - 1)*t22)
          where (mpg%cvNv(mpg%cvNvP(iCv2,1):
     &                    mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) == iCv)
           aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) =
     &      aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1)-
     &      t3*t2
          end where
         enddo
        endif
      enddo
*   ..check for empty matrix lines and replace with dummy equation
! $AD II-LOOP
      do iCv = 1, nCv
        lnonzero=.false.
        do inCv = 1,mpg%cvNvP(iCv,2)
          lnonzero = lnonzero.or.(aa(mpg%cvNvP(iCv,1)+inCv-1).ne.0.0_R8)
        enddo
        if (.not.lnonzero) then
          aa(mpg%cvNvP(iCv,1)) = 1.0_R8
          reshi(iCv) = 0.0_R8
        endif
      enddo
      if (switch%b2usht_iout.ne.0 .and. .not. in_parallel()) then
       call my_out_us(70,mpg%nCmxNv,0,aa,'b2usht_aa_ti')
      endif

*   ..solve the correction equation
      call b2uxus (nCv, mpg, aa, itcnt, reshi, corti, 'b2usht_Ti')
#if (defined(_OPENMP) && !defined(NO_OPENMP_B2USHT))
!svdk don't update last_solve_9, not yet supported
#else 
      last_solve_9(0:mpg%nnreg(0)) = solveei(0:mpg%nnreg(0))
! $AD II-LOOP
      do iCv = 1, nCv
        corti(iCv) = corti(iCv)
      enddo
#endif
      endif                                                          !srv 22.05.18

      if (switch%b2usht_iout.ne.0 .and. .not. in_parallel()) then
       call my_out_us(70,mpg%nCmxNv,0,aa,'b2usht_hi_aa')
       call my_out_us(70,nCv,0,reshi,'b2usht_reshi')
       call my_out_us(70,nCv,0,corti,'b2usht_corti')
      endif

#if (defined(_OPENMP) && !defined(NO_OPENMP_B2USHT))
!!$omp end task
!svdk task 4: neutral temperature equation
!!$omp task
!$omp section
#endif      
      if (switch%tn_style.eq.2) then
* ..neutral atom temperature correction equation
      if (any(solveen(0:mpg%nnreg(0)))) then                                !srv 22.05.18
#ifdef _OPENMP 
      call restart_ma28_for_us
#else
      new_matrix = .false.
      do ireg = 0, mpg%nnreg(0)
        new_matrix = new_matrix.or.
     &   solveen(ireg).neqv.solveei(ireg)
      end do
      if (new_matrix.or..true.) call restart_ma28_for_us
#endif
*   ..compute flo0, con0
      flo0        = flon
      con0(:,:,0) = conn
*   ..initialize matrix to zero
      aa = 0.0_R8
*   ..compute the nine-point matrix
! $AD II-LOOP
      do iCv = 1, nCv
        ireg = mpg%cvReg(iCv)
        if (solveen(ireg)) then
          aa(mpg%cvNvP(iCv,1)) =
     &    -shn(iCv,1)-shn(iCv,3)*nn(iCv)
          if (switch%b2mndt_style.lt.2) then
            aa(mpg%cvNvP(iCv,1)) = aa(mpg%cvNvP(iCv,1))
     &     +abs(reshn(iCv))/(switch%b2npht_rxg*tn(iCv))
          endif
        endif
      enddo
! $AD II-LOOP
      do iFc = 1, nFc
        iCv1 = mpg%fcCv(iFc,1) ! cell number of c1
        iCv2 = mpg%fcCv(iFc,2) ! cell number of c2
*     ..contributions from flo0(:,0), con0(:,0,:) ('5-point' contributions)
        t1 = - 0.5_R8*flo0(iFc,0) - con0(iFc,0,0)
        t2 =   0.5_R8*flo0(iFc,0) - con0(iFc,0,0)
        if (solveen(mpg%cvReg(iCv1))) then
         aa(mpg%cvNvP(iCv1,1)) = aa(mpg%cvNvP(iCv1,1)) - t1
         where (mpg%cvNv(mpg%cvNvP(iCv1,1):
     &                   mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) == iCv2)
          aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) =
     &     aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) +
     &     t2
         end where
        endif
        if (solveen(mpg%cvReg(iCv2))) then
         aa(mpg%cvNvP(iCv2,1)) = aa(mpg%cvNvP(iCv2,1)) - t2
         where (mpg%cvNv(mpg%cvNvP(iCv2,1):
     &                   mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) == iCv1)
          aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) =
     &     aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) +
     &     t1
         end where
        endif
*     ..contributions from flo0(:,1), con0(:,1,:) ('9-point' contributions)
        iVx1 = mpg%fcVx(iFc,1) ! number of v1
        iVx2 = mpg%fcVx(iFc,2) ! number of v2
        t1 = - 0.5_R8*flo0(iFc,1) - con0(iFc,1,0)
        t2 =   0.5_R8*flo0(iFc,1) - con0(iFc,1,0)
        t11 = 0.0_R8
        do inCv = 1, mpg%vxCvP(iVx1,2)
          t11 = t11 + 1/geo%vxVol(mpg%vxCvP(iVx1,1)+inCv-1)
        enddo
        t22 = 0.0_R8
        do inCv = 1, mpg%vxCvP(iVx2,2)
          t22 = t22 + 1/geo%vxVol(mpg%vxCvP(iVx2,1)+inCv-1)
        enddo
        ! contribution to equation of iCv1
        if (solveen(mpg%cvReg(iCv1))) then
         ! contribution from iVx1
         do inCv = 1, mpg%vxCvP(iVx1,2)
          iCv = mpg%vxCv(mpg%vxCvP(iVx1,1) + inCv - 1)
          t3  = 1/(geo%vxVol(mpg%vxCvP(iVx1,1) + inCv - 1)*t11)
          where (mpg%cvNv(mpg%cvNvP(iCv1,1):
     &                    mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) == iCv)
           aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) =
     &      aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1)-
     &      t3*t1
          end where
         enddo
         ! contribution from iVx2
         do inCv = 1, mpg%vxCvP(iVx2,2)
          iCv = mpg%vxCv(mpg%vxCvP(iVx2,1) + inCv - 1)
          t3  = 1/(geo%vxVol(mpg%vxCvP(iVx2,1) + inCv - 1)*t22)
          where (mpg%cvNv(mpg%cvNvP(iCv1,1):
     &                    mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) == iCv)
           aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) =
     &      aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1)+
     &      t3*t2
          end where
         enddo
        endif
        ! contribution to equation of iCv2
        if (solveen(mpg%cvReg(iCv2))) then
         ! contribution from iVx1
         do inCv = 1, mpg%vxCvP(iVx1,2)
          iCv = mpg%vxCv(mpg%vxCvP(iVx1,1) + inCv - 1)
          t3  = 1/(geo%vxVol(mpg%vxCvP(iVx1,1) + inCv - 1)*t11)
          where (mpg%cvNv(mpg%cvNvP(iCv2,1):
     &                    mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) == iCv)
           aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) =
     &      aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1)+
     &      t3*t1
          end where
         enddo
         ! contribution from iVx2
         do inCv = 1, mpg%vxCvP(iVx2,2)
          iCv = mpg%vxCv(mpg%vxCvP(iVx2,1) + inCv - 1)
          t3  = 1/(geo%vxVol(mpg%vxCvP(iVx2,1) + inCv - 1)*t22)
          where (mpg%cvNv(mpg%cvNvP(iCv2,1):
     &                    mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) == iCv)
           aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) =
     &      aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1)-
     &      t3*t2
          end where
         enddo
        endif
      enddo
*   ..check for empty matrix lines and replace with dummy equation
! $AD II-LOOP
      do iCv = 1, nCv
        lnonzero=.false.
        do inCv = 1,mpg%cvNvP(iCv,2)
          lnonzero = lnonzero.or.(aa(mpg%cvNvP(iCv,1)+inCv-1).ne.0.0_R8)
        enddo
        if (.not.lnonzero) then
          aa(mpg%cvNvP(iCv,1)) = 1.0_R8
          reshn(iCv) = 0.0_R8
        endif
      enddo

*   ..solve the correction equation
      call b2uxus (nCv, mpg, aa, itcnt, reshn, cortn, 'b2usht_Tn')
#if (defined(_OPENMP) && !defined(NO_OPENMP_B2USHT))
!svdk don't update last_solve_9, not yet supported
#else 
      last_solve_9(0:mpg%nnreg(0)) = solveen(0:mpg%nnreg(0))
! $AD II-LOOP
      do iCv = 1, nCv
        cortn(iCv) = cortn(iCv)
      enddo
#endif

      endif                                                          !srv 22.05.18
      endif

      if (switch%b2usht_iout.ne.0 .and. .not. in_parallel()) then
        call my_out_us(70,mpg%nCmxNv,0,aa,'b2usht_hn_aa')
        call my_out_us(70,nCv,0,reshi,'b2usht_reshn')
        call my_out_us(70,nCv,0,corti,'b2usht_cortn')
      endif

#if (defined(_OPENMP) && !defined(NO_OPENMP_B2USHT))
!!$omp end task
!svdk task 5: turbulent kinetic energy equation
!!$omp task
!$omp section
#endif 
* ..turbulent kinetic energy correction equation
      corkt = 0.0_R8
      if (switch%solve_keps.gt.0.and.any(solvekt(0:mpg%nnreg(0)))) then
       call restart_ma28_for_us
*   ..compute flo0, con0
       flo0        = flokt
       con0(:,:,0) = conkt
*   ..initialize matrix to zero
       aa = 0.0_R8
*   ..compute the nine-point matrix
! $AD II-LOOP
       do iCv = 1, nCv
        ireg = mpg%cvReg(iCv)
        if (solvekt(ireg)) then
         aa(mpg%cvNvP(iCv,1)) =
     &    -skt(iCv,1)-skt(iCv,3)*ni(iCv,1)
         if (switch%b2mndt_style.lt.2) then
          aa(mpg%cvNvP(iCv,1)) = aa(mpg%cvNvP(iCv,1)) +
     &    abs(reskt(iCv))/(switch%b2npht_rxg*(kt(iCv)+1.0e-5_R8*ev))
         endif
        endif
       enddo
! $AD II-LOOP
       do iFc = 1, nFc
        iCv1 = mpg%fcCv(iFc,1) ! cell number of c1
        iCv2 = mpg%fcCv(iFc,2) ! cell number of c2
*     ..contributions from flo0(:,0), con0(:,0,:) ('5-point' contributions)
        t1 = - 0.5_R8*flo0(iFc,0) - con0(iFc,0,0)
        t2 =   0.5_R8*flo0(iFc,0) - con0(iFc,0,0)
        if (solvekt(mpg%cvReg(iCv1))) then
          aa(mpg%cvNvP(iCv1,1)) = aa(mpg%cvNvP(iCv1,1)) - t1
          where (mpg%cvNv(mpg%cvNvP(iCv1,1):
     &                  mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) == iCv2)
           aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1)=
     &      aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1)+
     &      t2
          end where
        endif
        if (solvekt(mpg%cvReg(iCv2))) then
          aa(mpg%cvNvP(iCv2,1)) = aa(mpg%cvNvP(iCv2,1)) - t2
          where (mpg%cvNv(mpg%cvNvP(iCv2,1):
     &                  mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) == iCv1)
           aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1)=
     &      aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1)+
     &      t1
          end where
        endif
*     ..contributions from flo0(:,1), con0(:,1,:) ('9-point' contributions)
        iVx1 = mpg%fcVx(iFc,1) ! number of v1
        iVx2 = mpg%fcVx(iFc,2) ! number of v2
        t1 = - 0.5_R8*flo0(iFc,1) - con0(iFc,1,0)
        t2 =   0.5_R8*flo0(iFc,1) - con0(iFc,1,0)
        t11 = 0.0_R8
        do inCv = 1, mpg%vxCvP(iVx1,2)
          t11 = t11 + 1/geo%vxVol(mpg%vxCvP(iVx1,1)+inCv-1)
        enddo
        t22 = 0.0_R8
        do inCv = 1, mpg%vxCvP(iVx2,2)
          t22 = t22 + 1/geo%vxVol(mpg%vxCvP(iVx2,1)+inCv-1)
        enddo
        ! contribution to equation of iCv1
        if (solvekt(mpg%cvReg(iCv1))) then
         ! contribution from iVx1
         do inCv = 1, mpg%vxCvP(iVx1,2)
          iCv = mpg%vxCv(mpg%vxCvP(iVx1,1) + inCv - 1)
          t3  = 1/(geo%vxVol(mpg%vxCvP(iVx1,1) + inCv - 1)*t11)
          where (mpg%cvNv(mpg%cvNvP(iCv1,1):
     &                  mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) == iCv)
           aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) =
     &      aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1)-
     &      t3*t1
          end where
         enddo
         ! contribution from iVx2
         do inCv = 1, mpg%vxCvP(iVx2,2)
          iCv = mpg%vxCv(mpg%vxCvP(iVx2,1) + inCv - 1)
          t3  = 1/(geo%vxVol(mpg%vxCvP(iVx2,1) + inCv - 1)*t22)
          where (mpg%cvNv(mpg%cvNvP(iCv1,1):
     &                  mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) == iCv)
           aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) =
     &      aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1)+
     &      t3*t2
          end where
         enddo
        endif
        ! contribution to equation of iCv2
        if (solvekt(mpg%cvReg(iCv2))) then
         ! contribution from iVx1
         do inCv = 1, mpg%vxCvP(iVx1,2)
          iCv = mpg%vxCv(mpg%vxCvP(iVx1,1) + inCv - 1)
          t3  = 1/(geo%vxVol(mpg%vxCvP(iVx1,1) + inCv - 1)*t11)
          where (mpg%cvNv(mpg%cvNvP(iCv2,1):
     &                  mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) == iCv)
           aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) =
     &      aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1)+
     &      t3*t1
          end where
         enddo
         ! contribution from iVx2
         do inCv = 1, mpg%vxCvP(iVx2,2)
          iCv = mpg%vxCv(mpg%vxCvP(iVx2,1) + inCv - 1)
          t3  = 1/(geo%vxVol(mpg%vxCvP(iVx2,1) + inCv - 1)*t22)
          where (mpg%cvNv(mpg%cvNvP(iCv2,1):
     &                  mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) == iCv)
           aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) =
     &      aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1)-
     &      t3*t2
          end where
         enddo
        endif
       enddo
*   ..check for empty matrix lines and replace with dummy equation
! $AD II-LOOP
       do iCv = 1, nCv
        lnonzero=.false.
        do inCv = 1,mpg%cvNvP(iCv,2)
          lnonzero = lnonzero.or.(aa(mpg%cvNvP(iCv,1)+inCv-1).ne.0.0_R8)
        enddo
        if (.not.lnonzero) then
          aa(mpg%cvNvP(iCv,1)) = 1.0_R8
          reskt(iCv) = 0.0_R8
        endif
       enddo

*   ..solve the correction equation
       call b2uxus (nCv, mpg, aa, itcnt, reskt, corkt, 'b2usht_kt')
! $AD II-LOOP
       do iCv = 1, nCv
        corkt(iCv) = corkt(iCv)
       enddo
      endif !solve_keps - kt

      if (switch%b2usht_iout.ne.0 .and. .not. in_parallel()) then
       call my_out_us(70,mpg%nCmxNv,0,aa,'b2usht_kt_aa')
       call my_out_us(70,nCv,0,reskt,'b2usht_reskt')
       call my_out_us(70,nCv,0,corkt,'b2usht_corkt')
      endif

#if (defined(_OPENMP) && !defined(NO_OPENMP_B2USHT))
!!$omp end task
!svdk task 6: turbulent enstrophy equation
!!$omp task
!$omp section
#endif 
* ..turbulent enstrophy correction equation
      corzt = 0.0_R8
      if (switch%solve_keps.gt.1.and.any(solvezt(0:mpg%nnreg(0)))) then
       call restart_ma28_for_us
*   ..compute flo0, con0
       flo0        = flozt
       con0(:,:,0) = conzt
*   ..initialize matrix to zero
       aa = 0.0_R8
*   ..compute the nine-point matrix
! $AD II-LOOP
       do iCv = 1, nCv
        ireg = mpg%cvReg(iCv)
        if (solvezt(ireg)) then
         aa(mpg%cvNvP(iCv,1)) =
     &    -szt(iCv,1)-szt(iCv,3)*ni(iCv,1)
         if (switch%b2mndt_style.lt.2) then
           aa(mpg%cvNvP(iCv,1)) = aa(mpg%cvNvP(iCv,1)) +
     &     abs(reszt(iCv))/(switch%b2npht_rxg*(zt(iCv)+1.0e-5_R8*ev))
         endif
        endif
       enddo
! $AD II-LOOP
       do iFc = 1, nFc
        iCv1 = mpg%fcCv(iFc,1) ! cell number of c1
        iCv2 = mpg%fcCv(iFc,2) ! cell number of c2
*     ..contributions from flo0(:,0), con0(:,0,:) ('5-point' contributions)
        t1 = - 0.5_R8*flo0(iFc,0) - con0(iFc,0,0)
        t2 =   0.5_R8*flo0(iFc,0) - con0(iFc,0,0)
        if (solvezt(mpg%cvReg(iCv1))) then
          aa(mpg%cvNvP(iCv1,1)) = aa(mpg%cvNvP(iCv1,1)) - t1
          where (mpg%cvNv(mpg%cvNvP(iCv1,1):
     &                  mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) == iCv2)
           aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1)=
     &      aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1)+
     &      t2
          end where
        endif
        if (solvezt(mpg%cvReg(iCv2))) then
          aa(mpg%cvNvP(iCv2,1)) = aa(mpg%cvNvP(iCv2,1)) - t2
          where (mpg%cvNv(mpg%cvNvP(iCv2,1):
     &                  mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) == iCv1)
           aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1)=
     &      aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1)+
     &      t1
          end where
        endif
*     ..contributions from flo0(:,1), con0(:,1,:) ('9-point' contributions)
        iVx1 = mpg%fcVx(iFc,1) ! number of v1
        iVx2 = mpg%fcVx(iFc,2) ! number of v2
        t1 = - 0.5_R8*flo0(iFc,1) - con0(iFc,1,0)
        t2 =   0.5_R8*flo0(iFc,1) - con0(iFc,1,0)
        t11 = 0.0_R8
        do inCv = 1, mpg%vxCvP(iVx1,2)
          t11 = t11 + 1/geo%vxVol(mpg%vxCvP(iVx1,1)+inCv-1)
        enddo
        t22 = 0.0_R8
        do inCv = 1, mpg%vxCvP(iVx2,2)
          t22 = t22 + 1/geo%vxVol(mpg%vxCvP(iVx2,1)+inCv-1)
        enddo
        ! contribution to equation of iCv1
        if (solvezt(mpg%cvReg(iCv1))) then
         ! contribution from iVx1
         do inCv = 1, mpg%vxCvP(iVx1,2)
          iCv = mpg%vxCv(mpg%vxCvP(iVx1,1) + inCv - 1)
          t3  = 1/(geo%vxVol(mpg%vxCvP(iVx1,1) + inCv - 1)*t11)
          where (mpg%cvNv(mpg%cvNvP(iCv1,1):
     &                  mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) == iCv)
           aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) =
     &      aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1)-
     &      t3*t1
          end where
         enddo
         ! contribution from iVx2
         do inCv = 1, mpg%vxCvP(iVx2,2)
          iCv = mpg%vxCv(mpg%vxCvP(iVx2,1) + inCv - 1)
          t3  = 1/(geo%vxVol(mpg%vxCvP(iVx2,1) + inCv - 1)*t22)
          where (mpg%cvNv(mpg%cvNvP(iCv1,1):
     &                  mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) == iCv)
           aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) =
     &      aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1)+
     &      t3*t2
          end where
         enddo
        endif
        ! contribution to equation of iCv2
        if (solvezt(mpg%cvReg(iCv2))) then
         ! contribution from iVx1
         do inCv = 1, mpg%vxCvP(iVx1,2)
          iCv = mpg%vxCv(mpg%vxCvP(iVx1,1) + inCv - 1)
          t3  = 1/(geo%vxVol(mpg%vxCvP(iVx1,1) + inCv - 1)*t11)
          where (mpg%cvNv(mpg%cvNvP(iCv2,1):
     &                  mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) == iCv)
           aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) =
     &      aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1)+
     &      t3*t1
          end where
         enddo
         ! contribution from iVx2
         do inCv = 1, mpg%vxCvP(iVx2,2)
          iCv = mpg%vxCv(mpg%vxCvP(iVx2,1) + inCv - 1)
          t3  = 1/(geo%vxVol(mpg%vxCvP(iVx2,1) + inCv - 1)*t22)
          where (mpg%cvNv(mpg%cvNvP(iCv2,1):
     &                  mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) == iCv)
           aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) =
     &      aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1)-
     &      t3*t2
          end where
         enddo
        endif
       enddo
*   ..check for empty matrix lines and replace with dummy equation
! $AD II-LOOP
       do iCv = 1, nCv
        lnonzero=.false.
        do inCv = 1,mpg%cvNvP(iCv,2)
          lnonzero = lnonzero.or.(aa(mpg%cvNvP(iCv,1)+inCv-1).ne.0.0_R8)
        enddo
        if (.not.lnonzero) then
          aa(mpg%cvNvP(iCv,1)) = 1.0_R8
          reszt(iCv) = 0.0_R8
        endif
       enddo

*   ..solve the correction equation
       call b2uxus (nCv, mpg, aa, itcnt, reszt, corzt, 'b2usht_zt')
! $AD II-LOOP
       do iCv = 1, nCv
        corzt(iCv) = corzt(iCv)
       enddo
      endif !solve_keps - enstrophy

      if (switch%b2usht_iout.ne.0 .and. .not. in_parallel()) then
       call my_out_us(70,mpg%nCmxNv,0,aa,'b2usht_zt_aa')
       call my_out_us(70,nCv,0,reszt,'b2usht_reszt')
       call my_out_us(70,nCv,0,corzt,'b2usht_corzt')
      endif
#if (defined(_OPENMP) && !defined(NO_OPENMP_B2USHT))
!$omp end sections
!!$omp end task
!!$omp end single
!$omp end parallel
#endif 

*   ..return
      ncall_b2usht = ncall_b2usht+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2usht

      end subroutine b2usht
