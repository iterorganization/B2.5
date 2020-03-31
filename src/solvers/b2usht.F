*-----------------------------------------------------------------------
*.specification

      subroutine b2usht (nCv, nFc, nVx, ns, geo, mpg, nregionv, itcnt,
     &  rxg, ne, ni, te, ti, tn, 
     &  floe, cone, floi, coni, flon, conn,
     &  pccm1, she, shi, shn, ceqp,
     &  resht, reshe, reshi, reshn,
     &  cortt, corte, corti, cortn, flo0, con0, aa)
      use b2mod_types
      use b2mod_numerics_namelist
      use b2us_geo
      use b2us_map
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, nregionv, itcnt
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *  rxg, ne(nCv), ni(nCv,0:1), te(nCv), ti(nCv), tn(nCv), 
     &  floe(nFc,0:1), floi(nFc,0:1), flon(nFc,0:1),
     &  cone(nFc,0:1,0:2), coni(nFc,0:1,0:2), conn(nFc,0:1,0:2),       !srv 02.01.07
     &  pccm1(nCv)
      real (kind=R8) ::
     *  she(nCv,0:3), shi(nCv,0:3), shn(nCv,0:3), ceqp(nCv), 
     &  resht(nCv), reshe(nCv), reshi(nCv), reshn(nCv)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   cortt(nCv), corte(nCv), corti(nCv), cortn(nCv)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *  flo0(nFc,0:1), con0(nFc,0:1,0:2), aa(mpg%nCmxNv)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     b2usht computes a correction to the fields te and ti.
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
      integer ncall, iCv, iFc, iCv1, iCv2, iVx1, iVx2, inCv
      integer, save :: no_Ptncr_x_he = 0, no_Ptncr_y_he = 0            !srv 06.07.06
      integer, save :: no_Ptncr_x_hi = 0, no_Ptncr_y_hi = 0            !srv 06.07.06
      integer, save :: iadd_9 = 1                                      !srv 02.01.07
      integer, save :: last_ion
      integer, save :: iout = 0
      integer, save :: tn_style = 0
      integer ireg
      logical lnonzero
      save ncall
      real (kind=R8) ::
     *   t0, t1, t2, t3, t11, t22
      logical new_matrix
*   ..procedures
      intrinsic abs, sqrt, dim
      real (kind=R8) ::
     *   scalr
      external subini, subend, xertst, sfill
      external b2xvsg, b2xvff, b2xvfx, b2uxus, b2uxm9                 !srv 02.01.07
      scalr(t0) = 1.0e0_R8  ! t0                                      !srv 02.01.07
*   ..initialization
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2usht')
      if (ncall.eq.0) then
       call ipgeti ('b2tfhe_no_Ptncr_x', no_Ptncr_x_he)               !srv 06.07.06
       call ipgeti ('b2tfhi_no_Ptncr_x', no_Ptncr_x_hi)               !srv 06.07.06
       call ipgeti ('b2tfhe_no_Ptncr_y', no_Ptncr_y_he)               !srv 06.07.06
       call ipgeti ('b2tfhi_no_Ptncr_y', no_Ptncr_y_hi)               !srv 06.07.06
       call ipgeti ('b2usht_add_9', iadd_9)                           !srv 02.01.07
       call ipgeti ('b2usht_iout', iout)
       call ipgeti ('b2mn_tn_style', tn_style)
      endif
*   ..test nCv, nFc, ns
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test itcnt, rxg
      call xertst (0.le.itcnt, 'faulty argument itcnt')
      call xertst (0.lt.rxg, 'faulty argument rxg')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of ne, ni, te, ti
       call b2xvsg (nCv, ne, 1, 'ne', '.gt.')
       call b2xvsg (nCv, ni, 1, 'ni', '.gt.')
       call b2xvsg (nCv, te, 1, 'te', '.gt.')
       call b2xvsg (nCv, ti, 1, 'ti', '.gt.')
*    ..test sign of she, shi, ceqp
       call b2xvsg (nCv, she(1,0), 1, 'she0', '.ge.')
       call b2xvsg (nCv, she(1,1), 1, 'she1', '.le.')
       call b2xvsg (nCv, she(1,2), 1, 'she2', '.ge.')
       call b2xvsg (nCv, she(1,3), 1, 'she3', '.le.')
       call b2xvsg (nCv, shi(1,0), 1, 'shi0', '.ge.')
       call b2xvsg (nCv, shi(1,1), 1, 'shi1', '.le.')
       call b2xvsg (nCv, shi(1,2), 1, 'shi2', '.ge.')
       call b2xvsg (nCv, shi(1,3), 1, 'shi3', '.le.')
       call b2xvsg (nCv, shn(1,0), 1, 'shn0', '.ge.')
       call b2xvsg (nCv, shn(1,1), 1, 'shn1', '.le.')
       call b2xvsg (nCv, shn(1,2), 1, 'shn2', '.ge.')
       call b2xvsg (nCv, shn(1,3), 1, 'shn3', '.le.')
       call b2xvsg (nCv, ceqp, 1, 'ceqp', '.ge.')
      endif

      if (tn_style.le.1) then !! in case we solve a combined ion-neutral energy equation (default) or pure ion energy equation we make sure
         conn = 0.0_R8       !! the neutral fluxes and sources are zero.
         flon = 0.0_R8
         shn = 0.0_R8
      endif

      if (any(solveee(0:mpg%nnreg(0))).and.
     &    any(solveei(0:mpg%nnreg(0))).and.
     &    any(solveet(0:mpg%nnreg(0)))) then                                 !srv 22.05.18
      new_matrix = .false.
      do ireg = 0, nregionv
        new_matrix = new_matrix.or.
     &   solveet(ireg).neqv.last_solve_9(ireg)
      end do
      if (new_matrix.or..true.) call restart_ma28_for_us
* ..total heat correction equation  !! AFN_TODO switches for when to use solveen. Make compatible with branch stefano
*   ..compute flo0, con0
      flo0        = floe + floi + flon
      con0(:,:,0) = cone(:,:,0) + coni(:,:,0)+ conn(:,:,0)
*   ..initialize matrix to zero
      aa = 0.0e0_R8
*   ..compute the nine-point matrix
      do iCv = 1, nCv
        ireg = mpg%cvReg(iCv)
        if (solveee(ireg) .and. solveei(ireg). and. solveet(ireg)) then !srv 22.05.18
          aa(mpg%cvNvP(iCv,1)) =
     &    -she(iCv,1)-she(iCv,3)*ne(iCv)-
     &     shi(iCv,1)-shi(iCv,3)*ni(iCv,0)-
     &     shn(iCv,1)-shn(iCv,3)*ni(iCv,0)
     &    +abs(resht(iCv))/(rxg*min(te(iCv),ti(iCv)))           !srv 02.01.07
        else if (solveee(ireg) .and. .not.solveei(ireg)) then
          aa(mpg%cvNvP(iCv,1)) =
     &    -she(iCv,1)-she(iCv,3)*ne(iCv) !! AFN_TODO
     &    +abs(resht(iCv))/(rxg*te(iCv))                          !srv 02.01.07
        else if (.not.solveee(ireg) .and. solveei(ireg)) then
          aa(mpg%cvNvP(iCv,1)) =
     &    -shi(iCv,1)-shi(iCv,3)*ni(iCv,0) !! AFN_TODO
     &    +abs(resht(iCv))/(rxg*ti(iCv))                          !srv 02.01.07
        endif
      enddo
      do iFc = 1, nFc
        iCv1 = mpg%fcCv(iFc,1) ! cell number of c1
        iCv2 = mpg%fcCv(iFc,2) ! cell number of c2
*     ..contributions from flo0(:,0), con0(:,0,:) ('5-point' contributions)
        t1 = - 0.5_R8*flo0(iFc,0) - con0(iFc,0,0)
        t2 =   0.5_R8*flo0(iFc,0) - con0(iFc,0,0) 
        if (solveee(mpg%cvReg(iCv1)).or.solveei(mpg%cvReg(iCv1))) then
         aa(mpg%cvNvP(iCv1,1)) = aa(mpg%cvNvP(iCv1,1)) - t1
         where (mpg%cvNv(mpg%cvNvP(iCv1,1):
     &                   mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) == iCv2)
          aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) = 
     &     aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) +
     &     t2
         end where
        endif
        if (solveee(mpg%cvReg(iCv2)).or.solveei(mpg%cvReg(iCv2))) then
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
        if (solveee(mpg%cvReg(iCv1)).or.solveei(mpg%cvReg(iCv1))) then
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
        if (solveee(mpg%cvReg(iCv2)).or.solveei(mpg%cvReg(iCv2))) then
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
*   ..add matrix from contribution nine point stencil
#ifdef WG_TODO
      if(iadd_9 .eq. 1) then
       call b2uxm9 (nx, ny, nregionv, itcnt, solveee, cone(-1,-1,0,0,1),
     &  aa, 'b2usht_totale')                                         !srv 02.01.07
      endif
      if(iadd_9 .eq. 1) then
       call b2uxm9 (nx, ny, nregionv, itcnt, solveei, coni(-1,-1,0,0,1),
     &  aa, 'b2usht_totali')                                         !srv 02.01.07
      endif
#endif
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

*   ..solve the correction equation
      call b2uxus (nCv, mpg, aa, itcnt, resht, cortt, 'b2usht_total')  !srv 18.05.02
      endif                                                              !srv 22.05.18
      
      if (iout.ne.0) then
#ifdef WG_TODO
       call my_out_us(70,nCv,0,aa(-1,-1,-6),'b2usht_total_aam6')
       call my_out_us(70,nCv,0,aa(-1,-1,-5),'b2usht_total_aam5')
       call my_out_us(70,nCv,0,aa(-1,-1,-4),'b2usht_total_aam4')
       call my_out_us(70,nCv,0,aa(-1,-1,-3),'b2usht_total_aam3')
       call my_out_us(70,nCv,0,aa(-1,-1,-2),'b2usht_total_aam2')
       call my_out_us(70,nCv,0,aa(-1,-1,-1),'b2usht_total_aam1')
       call my_out_us(70,nCv,0,aa(-1,-1, 0),'b2usht_total_aam0')
       call my_out_us(70,nCv,0,aa(-1,-1, 1),'b2usht_total_aap1')
       call my_out_us(70,nCv,0,aa(-1,-1, 2),'b2usht_total_aap2')
       call my_out_us(70,nCv,0,aa(-1,-1, 3),'b2usht_total_aap3')
       call my_out_us(70,nCv,0,aa(-1,-1, 4),'b2usht_total_aap4')
       call my_out_us(70,nCv,0,aa(-1,-1, 5),'b2usht_total_aap5')
       call my_out_us(70,nCv,0,aa(-1,-1, 6),'b2usht_total_aap6')
#endif
       call my_out_us(70,nCv,0,resht,'b2usht_resht')
       call my_out_us(70,nCv,0,cortt,'b2usht_cortt')
      endif

* ..electron temperature correction equation
      if (any(solveee(0:mpg%nnreg(0)))) then                                 !srv 22.05.18
      new_matrix = .false.
      do ireg = 0, nregionv
        new_matrix = new_matrix.or.
     &   solveee(ireg).neqv.solveet(ireg)
      end do
      if (new_matrix.or..true.) call restart_ma28_for_us
*   ..compute flo0, con0
      flo0        = floe
      con0(:,:,0) = cone(:,:,0)
*   ..initialize matrix to zero
      aa = 0.0e0_R8
*   ..compute the nine-point matrix
      do iCv = 1, nCv
        ireg = mpg%cvReg(iCv)
        if (solveee(ireg)) then
          aa(mpg%cvNvP(iCv,1)) =
     &    -she(iCv,1)-she(iCv,3)*ne(iCv)+
     &    dim(1.0_R8,1.5_R8*(te(iCv)-ti(iCv))/te(iCv))*               !srv 02.01.07
     &    ceqp(iCv)*ne(iCv)*ni(iCv,0)/(te(iCv)*sqrt(te(iCv)))
     &    +abs(reshe(iCv))/(rxg*te(iCv))                              !srv 02.01.07
        endif
      enddo
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
*   ..add matrix from contribution nine point stencil
#ifdef WG_TODO
      if(iadd_9 .eq. 1) then
       call b2uxm9 (nx, ny, nregionv, itcnt, solveee, cone(-1,-1,0,0,1),
     &  aa, 'b2usht_Te')                                             !srv 02.01.07
      endif
#endif
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

*   ..solve the correction equation
      call b2uxus (nCv, mpg, aa, itcnt, reshe, corte, 'b2usht_Te')     !srv 18.05.02
      endif                                                              !srv 22.05.18

      if (iout.ne.0) then
#ifdef WG_TODO
       call my_out_us(70,nCv,0,aa(-1,-1,-6),'b2usht_he_aam6')
       call my_out_us(70,nCv,0,aa(-1,-1,-5),'b2usht_he_aam5')
       call my_out_us(70,nCv,0,aa(-1,-1,-4),'b2usht_he_aam4')
       call my_out_us(70,nCv,0,aa(-1,-1,-3),'b2usht_he_aam3')
       call my_out_us(70,nCv,0,aa(-1,-1,-2),'b2usht_he_aam2')
       call my_out_us(70,nCv,0,aa(-1,-1,-1),'b2usht_he_aam1')
       call my_out_us(70,nCv,0,aa(-1,-1, 0),'b2usht_he_aam0')
       call my_out_us(70,nCv,0,aa(-1,-1, 1),'b2usht_he_aap1')
       call my_out_us(70,nCv,0,aa(-1,-1, 2),'b2usht_he_aap2')
       call my_out_us(70,nCv,0,aa(-1,-1, 3),'b2usht_he_aap3')
       call my_out_us(70,nCv,0,aa(-1,-1, 4),'b2usht_he_aap4')
       call my_out_us(70,nCv,0,aa(-1,-1, 5),'b2usht_he_aap5')
       call my_out_us(70,nCv,0,aa(-1,-1, 6),'b2usht_he_aap6')
#endif
       call my_out_us(70,nCv,0,reshe,'b2usht_reshe')
       call my_out_us(70,nCv,0,corte,'b2usht_corte')
      endif

* ..atom temperature correction equation
      if (any(solveei(0:mpg%nnreg(0)))) then                                 !srv 22.05.18
      new_matrix = .false.
      do ireg = 0, nregionv
        new_matrix = new_matrix.or.
     &   solveei(ireg).neqv.solveee(ireg)
      end do
      if (new_matrix.or..true.) call restart_ma28_for_us
*   ..compute flo0, con0
      flo0        = floi
      con0(:,:,0) = coni(:,:,0)
*   ..initialize matrix to zero
      aa = 0.0e0_R8
*   ..compute the nine-point matrix
      do iCv = 1, nCv
        ireg = mpg%cvReg(iCv)
        if (solveei(ireg)) then
          aa(mpg%cvNvP(iCv,1)) =
     &    (-shi(iCv,1)-shi(iCv,3)*ni(iCv,0)+
     &    ceqp(iCv)*ne(iCv)*ni(iCv,0)/(te(iCv)*sqrt(te(iCv)))
     &    +abs(reshi(iCv))/(rxg*ti(iCv)))*scalr(ti(iCv))                              !srv 02.01.07
        endif
      enddo
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
*   ..add matrix from contribution nine point stencil
#ifdef WG_TODO
      if(iadd_9 .eq. 1) then
       call b2uxm9 (nx, ny, nregionv, itcnt, solveei, coni(-1,-1,0,0,1),
     &  aa, 'b2usht_Ti')                                             !srv 02.01.07
      endif
#endif
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

*   ..solve the correction equation
      call b2uxus (nCv, mpg, aa, itcnt, reshi, corti, 'b2usht_Ti')
      last_solve_9(0:mpg%nnreg(0)) = solveei(0:mpg%nnreg(0))
      do iCv = 1, nCv
        corti(iCv) = corti(iCv)*scalr(ti(iCv))
      enddo
      endif                                                          !srv 22.05.18


      if (iout.ne.0) then
#ifdef WG_TODO
       call my_out_us(70,nCv,0,aa(-1,-1,-6),'b2usht_hi_aam6')
       call my_out_us(70,nCv,0,aa(-1,-1,-5),'b2usht_hi_aam5')
       call my_out_us(70,nCv,0,aa(-1,-1,-4),'b2usht_hi_aam4')
       call my_out_us(70,nCv,0,aa(-1,-1,-3),'b2usht_hi_aam3')
       call my_out_us(70,nCv,0,aa(-1,-1,-2),'b2usht_hi_aam2')
       call my_out_us(70,nCv,0,aa(-1,-1,-1),'b2usht_hi_aam1')
       call my_out_us(70,nCv,0,aa(-1,-1, 0),'b2usht_hi_aam0')
       call my_out_us(70,nCv,0,aa(-1,-1, 1),'b2usht_hi_aap1')
       call my_out_us(70,nCv,0,aa(-1,-1, 2),'b2usht_hi_aap2')
       call my_out_us(70,nCv,0,aa(-1,-1, 3),'b2usht_hi_aap3')
       call my_out_us(70,nCv,0,aa(-1,-1, 4),'b2usht_hi_aap4')
       call my_out_us(70,nCv,0,aa(-1,-1, 5),'b2usht_hi_aap5')
       call my_out_us(70,nCv,0,aa(-1,-1, 6),'b2usht_hi_aap6')
#endif
       call my_out_us(70,nCv,0,reshi,'b2usht_reshi')
       call my_out_us(70,nCv,0,corti,'b2usht_corti')
      endif

* ..neutral atom temperature correction equation
      if (any(solveen(0:mpg%nnreg(0)))) then                                 !srv 22.05.18
      new_matrix = .false.
      do ireg = 0, nregionv
        new_matrix = new_matrix.or.
     &   solveen(ireg).neqv.solveei(ireg)
      end do
      if (new_matrix.or..true.) call restart_ma28_for_us
*   ..compute flo0, con0
      flo0        = flon
      con0(:,:,0) = conn(:,:,0)
*   ..initialize matrix to zero
      aa = 0.0e0_R8
*   ..compute the nine-point matrix
      do iCv = 1, nCv
        ireg = mpg%cvReg(iCv)
        if (solveen(ireg)) then
          aa(mpg%cvNvP(iCv,1)) =
     &    (-shn(iCv,1)-shn(iCv,3)*ni(iCv,0)+
     &    +abs(reshn(iCv))/(rxg*tn(iCv)))*scalr(tn(iCv))                              !srv 02.01.07
        endif
      enddo
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
*   ..add matrix from contribution nine point stencil
#ifdef WG_TODO
      if(iadd_9 .eq. 1) then
       call b2uxm9 (nx, ny, nregionv, itcnt, solveen, conn(-1,-1,0,0,1),
     &  aa, 'b2usht_Tn')                                             !srv 02.01.07
      endif
#endif
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
      last_solve_9(0:mpg%nnreg(0)) = solveen(0:mpg%nnreg(0))
      do iCv = 1, nCv
        cortn(iCv) = cortn(iCv)*scalr(tn(iCv))
      enddo
      endif                                                          !srv 22.05.18


      if (iout.ne.0) then
#ifdef WG_TODO
       call my_out_us(70,nCv,0,aa(-1,-1,-6),'b2usht_hn_aam6')
       call my_out_us(70,nCv,0,aa(-1,-1,-5),'b2usht_hn_aam5')
       call my_out_us(70,nCv,0,aa(-1,-1,-4),'b2usht_hn_aam4')
       call my_out_us(70,nCv,0,aa(-1,-1,-3),'b2usht_hn_aam3')
       call my_out_us(70,nCv,0,aa(-1,-1,-2),'b2usht_hn_aam2')
       call my_out_us(70,nCv,0,aa(-1,-1,-1),'b2usht_hn_aam1')
       call my_out_us(70,nCv,0,aa(-1,-1, 0),'b2usht_hn_aam0')
       call my_out_us(70,nCv,0,aa(-1,-1, 1),'b2usht_hn_aap1')
       call my_out_us(70,nCv,0,aa(-1,-1, 2),'b2usht_hn_aap2')
       call my_out_us(70,nCv,0,aa(-1,-1, 3),'b2usht_hn_aap3')
       call my_out_us(70,nCv,0,aa(-1,-1, 4),'b2usht_hn_aap4')
       call my_out_us(70,nCv,0,aa(-1,-1, 5),'b2usht_hn_aap5')
       call my_out_us(70,nCv,0,aa(-1,-1, 6),'b2usht_hn_aap6')
#endif
       call my_out_us(70,nCv,0,reshi,'b2usht_reshn')
       call my_out_us(70,nCv,0,corti,'b2usht_cortn')
      endif



*   ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2usht

      end
