*-----------------------------------------------------------------------
*.specification

      subroutine b2uspo (nCv, nFc, nVx, nregionv, solvereg, itcnt, 
     &  switch, geo,
     &  mpg, ne, te, conc, sch, respo, corpo, aa, name)                  !srv 18.03.10
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2uspo
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, nregionv, itcnt
      character*(*) name                                             !srv 18.03.10
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *  ne(nCv), te(nCv), conc(nFc,0:1), sch(nCv,0:3), respo(nCv)
      logical solvereg(0:nregionv)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  corpo(nCv)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *  aa(mpg%nCmxNv)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2USPO computes a correction towards the solution of the
*     potential equation.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, iFc, iCv1, iCv2, iVx1, iVx2, inCv, ireg, k
      character chark*1
      real (kind=R8) :: t1, t2, t3, t11, t22
      logical lnonzero
*   ..procedures
      external xertst, sfill
      external b2xvsg, b2xvff, b2ux5p

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2uspo')
*   ..test nCv, nFc
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
*   ..test itcnt
      call xertst (0.le.itcnt, 'faulty argument itcnt')
*   ..extensive tests on first few calls
      if (ncall_b2uspo.lt.3) then
*    ..test sign of ne, te, conc, sch
       call b2xvsg (nCv, ne, 1, 'ne', '.gt.')
       call b2xvsg (nCv, te, 1, 'te', '.gt.')
       call b2xvsg (nCv, sch(1,1), 1, 'sch1', '.le.')
       call b2xvsg (nCv, sch(1,3), 1, 'sch3', '.le.')
      endif

* ..compute the correction
      aa = 0.0_R8
*   ..compute the nine-point matrix
! $AD II-LOOP
      do iCv = 1, nCv
        ireg = mpg%cvReg(iCv)
        if (solvereg(ireg)) then
         aa(mpg%cvNvP(iCv,1)) = -sch(iCv,1)-sch(iCv,3)*ne(iCv)
        else
         aa(mpg%cvNvP(iCv,1)) = 1.0_R8
        endif
      enddo
! $AD II-LOOP
      do iFc = 1, nFc
        iCv1 = mpg%fcCv(iFc,1) ! cell number of c1
        iCv2 = mpg%fcCv(iFc,2) ! cell number of c2
*     ..contributions from conc(:,0) ('5-point' contributions)
        t1 = -conc(iFc,0)
        t2 = -conc(iFc,0) 
        if (solvereg(mpg%cvReg(iCv1))) then
         aa(mpg%cvNvP(iCv1,1)) = aa(mpg%cvNvP(iCv1,1)) - t1
         where (mpg%cvNv(mpg%cvNvP(iCv1,1):
     &                   mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) == iCv2)
          aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) = 
     &     aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) +
     &     t2
         end where
        endif
        if (solvereg(mpg%cvReg(iCv2))) then
         aa(mpg%cvNvP(iCv2,1)) = aa(mpg%cvNvP(iCv2,1)) - t2
         where (mpg%cvNv(mpg%cvNvP(iCv2,1):
     &                   mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) == iCv1)
          aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) = 
     &     aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) +
     &     t1 
         end where
        endif
*     ..contributions from conc(:,1) ('9-point' contributions)
        iVx1 = mpg%fcVx(iFc,1) ! number of v1
        iVx2 = mpg%fcVx(iFc,2) ! number of v2
        t1 = -conc(iFc,1)
        t2 = -conc(iFc,1)
        t11 = 0.0_R8
        do inCv = 1, mpg%vxCvP(iVx1,2)
          t11 = t11 + 1/geo%vxVol(mpg%vxCvP(iVx1,1)+inCv-1)
        enddo
        t22 = 0.0_R8
        do inCv = 1, mpg%vxCvP(iVx2,2)
          t22 = t22 + 1/geo%vxVol(mpg%vxCvP(iVx2,1)+inCv-1)
        enddo
        ! contribution to equation of iCv1
        if (solvereg(mpg%cvReg(iCv1))) then
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
        if (solvereg(mpg%cvReg(iCv2))) then
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
*   ..correct RHS if not solving
! $AD II-LOOP
      do iCv = 1, nCv
        lnonzero=.false.
        do inCv = 1,mpg%cvNvP(iCv,2)
          lnonzero = lnonzero.or.(aa(mpg%cvNvP(iCv,1)+inCv-1).ne.0.0_R8)
        enddo
        if (.not.lnonzero) then
          aa(mpg%cvNvP(iCv,1)) = 1.0_R8
          respo(iCv) = 0.0_R8
        endif
        if(.not.solvereg(mpg%cvReg(iCv))) respo(iCv) = 0.0_R8
      enddo

*   ..solve the correction equation
      write(*,'(1x,a,a)') 'Calling b2uxus from '//trim(name)
      call b2uxus (nCv, mpg, aa, itcnt, respo, corpo, name)           !srv 18.03.10

      if (switch%b2uspo_iout.ne.0) then                                          !srv 08.07.08 {
        do k=0,3
          write (chark,'(i1)') k
          call my_out_us(70,nCv,0,sch(1,k),'b2uspo_sch'//chark)
        enddo
#ifdef WG_TODO
        call my_out_us(70,nCv,0,aa(-1,-1,-6),'b2uspo_aam6')
        call my_out_us(70,nCv,0,aa(-1,-1,-5),'b2uspo_aam5')
        call my_out_us(70,nCv,0,aa(-1,-1,-4),'b2uspo_aam4')
        call my_out_us(70,nCv,0,aa(-1,-1,-3),'b2uspo_aam3')
        call my_out_us(70,nCv,0,aa(-1,-1,-2),'b2uspo_aam2')
        call my_out_us(70,nCv,0,aa(-1,-1,-1),'b2uspo_aam1')
        call my_out_us(70,nCv,0,aa(-1,-1, 0),'b2uspo_aam0')
        call my_out_us(70,nCv,0,aa(-1,-1, 1),'b2uspo_aap1')
        call my_out_us(70,nCv,0,aa(-1,-1, 2),'b2uspo_aap2')
        call my_out_us(70,nCv,0,aa(-1,-1, 3),'b2uspo_aap3')
        call my_out_us(70,nCv,0,aa(-1,-1, 4),'b2uspo_aap4')
        call my_out_us(70,nCv,0,aa(-1,-1, 5),'b2uspo_aap5')
        call my_out_us(70,nCv,0,aa(-1,-1, 6),'b2uspo_aap6')
#endif
        call my_out_us(70,nCv,0,respo,'b2uspo_respo')
        call my_out_us(70,nCv,0,corpo,'b2uspo_corpo')
      endif

* ..return
      ncall_b2uspo = ncall_b2uspo+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2uspo

      end subroutine b2uspo
