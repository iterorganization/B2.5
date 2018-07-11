*-----------------------------------------------------------------------
*.specification

      subroutine b2uspo (nCv, nFc, nVx, nregionv, solvereg, itcnt, geo,
     &  mpg, ne, te, conc, sch, respo, corpo, aa, name)                  !srv 18.03.10
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_map
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, nregionv, itcnt
      character*(*) name                                             !srv 18.03.10
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *  ne(nCv), te(nCv), conc(nFc,0:1), sch(nCv,0:3), respo(nCv)
      logical solvereg(0:nregionv)      
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  corpo(nCv)
*   ..workspace arguments (unspecified on entry and on exit)
!WG_TODO: better definition of size aa still needed
      real (kind=R8) ::
     *  aa(9*nCv)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     b2uspo computes a correction towards the solution of the
*     potential equation.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, iCv, iFc, iCv1, iCv2, iNv1, iNv2, ireg, k
      integer, save :: iout = 0
      character chark*1
      real (kind=R8) :: t1, t2
      logical lnonzero
      save ncall
*   ..procedures
      intrinsic min
      external subini, subend, xertst, sfill
      external b2xvsg, b2xvff, b2ux5p
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2uspo')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti ('b2uspo_iout', iout)
      endif  
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
*   ..test itcnt
      call xertst (0.le.itcnt, 'faulty argument itcnt')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of ne, te, conc, sch
       call b2xvsg (nCv, ne, 1, 'ne', '.gt.')
       call b2xvsg (nCv, te, 1, 'te', '.gt.')
       call b2xvsg (nCv, sch(1,1), 1, 'sch1', '.le.')
       call b2xvsg (nCv, sch(1,3), 1, 'sch3', '.le.')
      endif

* ..compute the correction
      aa = 0.0_R8
*   ..compute the nine-point matrix
      do iCv = 1, nCv
        ireg = mpg%cvReg(iCv)
        if (solvereg(ireg)) then
         aa(mpg%aaP(iCv,1)) = -sch(iCv,1)-sch(iCv,3)*ne(iCv)
        else
         aa(mpg%aaP(iCv,1)) = 1.0_R8
        endif
      enddo
      do iFc = 1, nFc
        iCv1 = mpg%fcCv(iFc,1) ! cell number of c1
        iCv2 = mpg%fcCv(iFc,2) ! cell number of c2
*     ..contributions from conc(:,0) ('5-point' contributions)
        iNv1 = mpg%fcNv(iFc,1) ! number of c2 as seen from c1 
        iNv2 = mpg%fcNv(iFc,2) ! number of c1 as seen from c2
        t1 = -conc(iFc,0)
        t2 = -conc(iFc,0) 
        if (solvereg(mpg%cvReg(iCv1))) then
          aa(mpg%aaP(iCv1,1))        = aa(mpg%aaP(iCv1,1))        - t1
          aa(mpg%aaP(iCv1,1)+iNv1-1) = aa(mpg%aaP(iCv1,1)+iNv1-1) + t2
        endif
        if (solvereg(mpg%cvReg(iCv2))) then
          aa(mpg%aaP(iCv2,1))        = aa(mpg%aaP(iCv2,1))        - t2
          aa(mpg%aaP(iCv2,1)+iNv2-1) = aa(mpg%aaP(iCv2,1)+iNv2-1) + t1
        endif
*     ..contributions from conc(:,1) ('9-point' contributions)
!WG_TODO
      enddo
*   ..correct RHS if not solving
      do iCv = 1, nCv
        lnonzero=.false.
        do iNv1 = 1,mpg%aaP(iCv,2)
          lnonzero = lnonzero.or.(aa(mpg%aaP(iCv,1)+iNv1-1).ne.0.0_R8)
        enddo
        if (.not.lnonzero) then
          aa(mpg%aaP(iCv,1)) = 1.0_R8
          respo(iCv) = 0.0_R8
        endif
        if(.not.solvereg(mpg%cvReg(iCv))) respo(iCv) = 0.0_R8
      enddo

*   ..solve the correction equation
      write(*,*) 'Calling b2uxus from '//trim(name)
      call b2uxus (nCv, mpg, aa, itcnt, respo, corpo, name)           !srv 18.03.10
#ifdef WG_TODO
      if (iout.ne.0) then                                          !srv 08.07.08 {
        do k=0,3
          write (chark,'(i1)') k
          call my_out_us(70,nCv,sch(-1,-1,k),'b2uspo_sch'//chark)
        enddo
        call my_out_us(70,nCv,aa(-1,-1,-6),'b2uspo_aam6')
        call my_out_us(70,nCv,aa(-1,-1,-5),'b2uspo_aam5')
        call my_out_us(70,nCv,aa(-1,-1,-4),'b2uspo_aam4')
        call my_out_us(70,nCv,aa(-1,-1,-3),'b2uspo_aam3')
        call my_out_us(70,nCv,aa(-1,-1,-2),'b2uspo_aam2')
        call my_out_us(70,nCv,aa(-1,-1,-1),'b2uspo_aam1')
        call my_out_us(70,nCv,aa(-1,-1, 0),'b2uspo_aam0')
        call my_out_us(70,nCv,aa(-1,-1, 1),'b2uspo_aap1')
        call my_out_us(70,nCv,aa(-1,-1, 2),'b2uspo_aap2')
        call my_out_us(70,nCv,aa(-1,-1, 3),'b2uspo_aap3')
        call my_out_us(70,nCv,aa(-1,-1, 4),'b2uspo_aap4')
        call my_out_us(70,nCv,aa(-1,-1, 5),'b2uspo_aap5')
        call my_out_us(70,nCv,aa(-1,-1, 6),'b2uspo_aap6')
        call my_out_us(70,nCv,respo,'b2uspo_respo')
        call my_out_us(70,nCv,corpo,'b2uspo_corpo')
      endif
#endif
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2uspo

      end
