*-----------------------------------------------------------------------
*.specification

      subroutine b2usht_mol (nCv, nFc, nVx, geo, mpg,
     &  itcnt, nm0, tm0, flo, con, shm, reshtm,
     &  cortm, flo0, con0, aa, switch)
      use b2mod_types
      use b2mod_numerics_namelist
      use b2mod_switches
      use b2mod_constants
      use b2us_geo
      use b2us_map
      use b2mod_ma28_for_us
     , , only : restart_ma28_for_us
      use b2mod_openmp
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, itcnt
      type (switches) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *  nm0(nCv),  tm0(nCv), flo(nFc,0:1), con(nFc,0:1,0:2)
      real (kind=R8) ::
     *  shm(nCv,0:3), reshtm(nCv)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  cortm(nCv) 
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *  flo0(nFc,0:1), con0(nFc,0:1,0:2), aa(mpg%nCmxNv)
      real (kind=R8) :: dumm(1)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2USHT computes a correction to the fields tm
*
*.declarations

*   ..local variables
      integer :: iCv, iFc, iCv1, iCv2, iVx1, iVx2, inCv
      integer, save :: ncall_b2usht_mol = 0
      integer :: ireg
      logical :: lnonzero
      real (kind=R8) ::
     *   t0, t1, t2, t3, t11, t22
      logical new_matrix
*   ..procedures
      intrinsic abs, sqrt, dim
      external xertst, sfill
      external b2xvsg, b2uxus                                         
*   ..initialization

*-----------------------------------------------------------------------
*.computation

* ..preliminaries

      iCv1 = 1 !csc needed for correct initialization in adjoint AD
      iCv2 = 1
*   ..subprogram start-up calls
      call subini ('b2usht_mol')
*   ..test nCv, nFc, nVx
      call xertst (0.lt.nCv.and.0.lt.nFc.and.0.lt.nVx,
     .  'faulty argument nCv, nFc, nVx')
*   ..test itcnt, rxg
      call xertst (0.le.itcnt, 'faulty argument itcnt')
      call xertst (0.0_R8.lt.switch%b2npht_rxg, 'faulty argument rxg')
*   ..extensive tests on first few calls
      if (ncall_b2usht_mol.lt.3) then
*    ..test sign of ne, ni, te, ti
       call b2xvsg (nCv, nm0, 1, 'nm', '.gt.')
       call b2xvsg (nCv, tm0, 1, 'tm', '.gt.')
       call b2xvsg (nCv, shm(1,0), 1, 'shm0', '.ge.')
       call b2xvsg (nCv, shm(1,1), 1, 'shm1', '.le.')
       call b2xvsg (nCv, shm(1,2), 1, 'shm2', '.ge.')
       call b2xvsg (nCv, shm(1,3), 1, 'shm3', '.le.')
      endif


* ..heat correction equation
      new_matrix = .false.
* ..total heat correction equation
*   ..compute flo0, con0
      flo0        = flo
      con0(:,:,0) = con(:,:,0) 


*   ..initialize matrix to zero
      aa = 0.0_R8
*   ..compute the nine-point matrix
      do iCv = 1, nCv
        ireg = mpg%cvReg(iCv)
        t0 = 1.0e60_R8
        aa(mpg%cvNvP(iCv,1)) = aa(mpg%cvNvP(iCv,1))
     &   -shm(iCv,1)-shm(iCv,3)*nm0(iCv)
        t0 = min(t0,tm0(iCv))
      enddo

      do iFc = 1, nFc
        iCv1 = mpg%fcCv(iFc,1) ! cell number of c1
        iCv2 = mpg%fcCv(iFc,2) ! cell number of c2
*     ..contributions from flo0(:,0), con0(:,0,:) ('5-point' contributions)
        t1 = - 0.5_R8*flo0(iFc,0) - con0(iFc,0,0)
        t2 =   0.5_R8*flo0(iFc,0) - con0(iFc,0,0)
        aa(mpg%cvNvP(iCv1,1)) = aa(mpg%cvNvP(iCv1,1)) - t1
        where (mpg%cvNv(mpg%cvNvP(iCv1,1):
     &                   mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) == iCv2)
          aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) =
     &     aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) +
     &     t2
        end where
        aa(mpg%cvNvP(iCv2,1)) = aa(mpg%cvNvP(iCv2,1)) - t2
        where (mpg%cvNv(mpg%cvNvP(iCv2,1):
     &                   mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) == iCv1)
          aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) =
     &     aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) +
     &     t1
        end where

*     ..contributions from flo0(:,1), con0(:,1,:) ('9-point' contributions)
        iVx1 = mpg%fcVx(iFc,1) ! number of v1
        iVx2 = mpg%fcVx(iFc,2) ! number of v2
        t1 = - 0.5_R8*flo0(iFc,1) - con0(iFc,1,0)
        t2 =   0.5_R8*flo0(iFc,1) - con0(iFc,1,0)
        t11 = 0.0_R8
        do inCv = 1, mpg%vxCvP(iVx1,2)
          t11 = t11 + 1.0_R8/geo%vxVol(mpg%vxCvP(iVx1,1)+inCv-1)
        enddo
        t22 = 0.0_R8
        do inCv = 1, mpg%vxCvP(iVx2,2)
          t22 = t22 + 1.0_R8/geo%vxVol(mpg%vxCvP(iVx2,1)+inCv-1)
        enddo
        ! contribution to equation of iCv1

        ! contribution from iVx1
        do inCv = 1, mpg%vxCvP(iVx1,2)
         iCv = mpg%vxCv(mpg%vxCvP(iVx1,1) + inCv - 1)
         t3  = 1.0_R8/(geo%vxVol(mpg%vxCvP(iVx1,1) + inCv - 1)*t11)
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
         t3  = 1.0_R8/(geo%vxVol(mpg%vxCvP(iVx2,1) + inCv - 1)*t22)
         where (mpg%cvNv(mpg%cvNvP(iCv1,1):
     &                    mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) == iCv)
           aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) =
     &      aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1)+
     &      t3*t2
         end where
        enddo

        ! contribution to equation of iCv2
        ! contribution from iVx1
        do inCv = 1, mpg%vxCvP(iVx1,2)
         iCv = mpg%vxCv(mpg%vxCvP(iVx1,1) + inCv - 1)
         t3  = 1.0_R8/(geo%vxVol(mpg%vxCvP(iVx1,1) + inCv - 1)*t11)
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
         t3  = 1.0_R8/(geo%vxVol(mpg%vxCvP(iVx2,1) + inCv - 1)*t22)
         where (mpg%cvNv(mpg%cvNvP(iCv2,1):
     &                    mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) == iCv)
           aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) =
     &      aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1)-
     &      t3*t2
         end where
        enddo
      enddo
*   ..check for empty matrix lines and replace with dummy equation
      do iCv = 1, nCv
        lnonzero=.false.
        do inCv = 1, mpg%cvNvP(iCv,2)
          lnonzero = lnonzero.or.(aa(mpg%cvNvP(iCv,1)+inCv-1).ne.0.0_R8)
        enddo
        if (.not.lnonzero) then
          aa(mpg%cvNvP(iCv,1)) = 1.0_R8
          reshtm(iCv) = 0.0_R8
        endif
      enddo
*   ..solve the correction equation
      call b2uxus (nCv, mpg, aa, itcnt, reshtm, cortm, 'b2usht_mol',  
     .             switch%b2uxus_style)

*   ..return
      ncall_b2usht_mol = ncall_b2usht_mol+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2usht_mol

      end subroutine b2usht_mol
