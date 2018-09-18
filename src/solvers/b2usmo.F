*-----------------------------------------------------------------------
*.specification

      subroutine b2usmo (nCv, nFc, nVx, isb, geo, mpg, 
     &  nregionv, solvereg, itcnt, rxg, rob, pb, ub,                !srv 09.01.01
     &  smb, flcb, cvsb, resmb, ctcfb, corub, pccb, aa, aad0, name) !srv 06.04.07
      use b2mod_types
      use b2us_geo
      use b2us_map
      implicit none
*   ..input arguments (unchaninvged on exit)
      integer nCv, nFc, nVx, isb, nregionv, itcnt
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      character*(*) name                                             !srv 18.05.02
      real (kind=R8) ::
     *  rxg, rob(nCv),
     &  pb(nCv), ub(nCv), smb(nCv,0:3), flcb(nFc,0:1), cvsb(nFc,0:1),
     &  resmb(nCv), ctcfb(nCv)                                       !srv 06.04.07
      logical solvereg(0:nregionv)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  corub(nCv), pccb(nCv,0:1)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *  aa(mpg%nCmxNv), aad0(nCv), flcb0(nFc,0:1), cvsb0(nFc,0:1)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     b2usmo computes a correction towards the solution of the
*     parallel momentum balance equation for one atomic species.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the velocity correction, corub, for one
*     species.
*     ...
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, iCv, iFc, iCv1, iCv2, iNv1, iNv2, ireg
      integer, save :: iout = 0                                       !srv 17.06.02
      real (kind=R8) ::
     *  cfc0, t0, t1, t2, t3, t4, flo, dif, flo9, con9, csb, csb1
      logical lnonzero
      save ncall, cfc0
*   ..procedures
      intrinsic max, abs, sqrt
      external subini, subend, xertst, ipgetr, sfill
      external b2xvsg, b2xvff, b2ux5p, intfaceh, intfacev
*   ..initialisation
      data ncall/0/, cfc0/1.0_R8/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2usmo')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2usmo_cfc0', cfc0)
       call xertst (0.0_R8.le.cfc0, 'faulty internal parameter cfc0')
       call ipgeti ('b2usmo_iout', iout)                              !srv 17.06.02
      endif
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
*   ..test itcnt
      call xertst (0.le.itcnt, 'faulty argument itcnt')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of rob, pb, smb1, smb3, cvsb
       call b2xvsg (nCv, rob, 1, 'rob', '.gt.')
       call b2xvsg (nCv, pb, 1, 'pb', '.gt.')
       call b2xvsg (nCv, smb(1,1), 1, 'smb1', '.le.')
       call b2xvsg (nCv, smb(1,3), 1, 'smb3', '.le.')
      endif

* ..apply discretization scheme
      call calccoef(nCv,nFc,nVx,2,geo,mpg,flcb,cvsb,flcb0,cvsb0)

* ..compute the correction
      aa = 0.0_R8
*   ..compute the nine-point matrix
      do iCv = 1, nCv
        ireg = mpg%cvReg(iCv)
        if (solvereg(ireg)) then
         csb = sqrt(pb(iCv)/rob(iCv))
         csb1 = sqrt(ub(iCv)**2+csb**2)
         aad0(iCv) = -smb(iCv,1)-
     &    smb(iCv,3)*rob(iCv)+
     &    cfc0*abs(geo%cvVol(iCv)*geo%cvH(iCv)/geo%cvHx(iCv)*       !srv 09.01.01
     &    geo%cvBb(iCv,0)/geo%cvBb(iCv,3)*csb)*rob(iCv)
         aa(mpg%aaP(iCv,1)) = aad0(iCv)+abs(resmb(iCv))/(rxg*csb1)
        else
         aa(mpg%aaP(iCv,1)) = 1.0_R8
         aad0(iCv) = 0.0_R8
        endif
      enddo
      do iFc = 1, nFc
        iCv1 = mpg%fcCv(iFc,1) ! cell number of c1
        iCv2 = mpg%fcCv(iFc,2) ! cell number of c2
*     ..contributions from flcb0(:,0), cvsb0(:,0) ('5-point' contributions)
        iNv1 = mpg%fcNv(iFc,1) ! number of c2 as seen from c1 
        iNv2 = mpg%fcNv(iFc,2) ! number of c1 as seen from c2
        t1 = - 0.5_R8*flcb0(iFc,0) - cvsb0(iFc,0)
        t2 =   0.5_R8*flcb0(iFc,0) - cvsb0(iFc,0)
        if (solvereg(mpg%cvReg(iCv1))) then
          aa(mpg%aaP(iCv1,1))        = aa(mpg%aaP(iCv1,1))        - t1
          aa(mpg%aaP(iCv1,1)+iNv1-1) = aa(mpg%aaP(iCv1,1)+iNv1-1) + t2
        endif
        if (solvereg(mpg%cvReg(iCv2))) then
          aa(mpg%aaP(iCv2,1))        = aa(mpg%aaP(iCv2,1))        - t2
          aa(mpg%aaP(iCv2,1)+iNv2-1) = aa(mpg%aaP(iCv2,1)+iNv2-1) + t1
        endif
*     ..contributions from flcb0(:,1), cvsb0(:,1) ('9-point' contributions)
!WG_TODO
      enddo
*   ..compute pccb
      pccb = 0.0_R8
      do iCv = 1, nCv
        ireg = mpg%cvReg(iCv)
        if (solvereg(ireg)) then
         pccb(iCv,0) = (geo%cvVol(iCv)*geo%cvH(iCv)/geo%cvHx(iCv))*      !srv 09.01.01
     .    geo%cvBb(iCv,0)/geo%cvBb(iCv,3)/aa(mpg%aaP(iCv,1))             !srv 11.01.13
         pccb(iCv,1) = (geo%cvVol(iCv)*geo%cvH(iCv)/geo%cvHx(iCv))*      !srv 09.01.01
     .    geo%cvBb(iCv,0)/geo%cvBb(iCv,3)/aad0(iCv)                      !srv 11.01.13
        else
         resmb(iCv) = 0.0_R8   ! take opportunity to correct RHS
        endif
      enddo
*   ..contribution from linearization of centrifugal force
      do iCv = 1, nCv
        aa(mpg%aaP(iCv,1)) = aa(mpg%aaP(iCv,1)) - 
     &    min(0.0e0_R8,ctcfb(iCv))
      enddo
      do iCv = 1, nCv
        lnonzero=.false.
        do iNv1 = 1,mpg%aaP(iCv,2)
          lnonzero = lnonzero.or.(aa(mpg%aaP(iCv,1)+iNv1-1).ne.0.0_R8)
        enddo
        if (.not.lnonzero) then
          aa(mpg%aaP(iCv,1)) = 1.0_R8
          resmb(iCv) = 0.0_R8
        endif
      enddo

*   ..solve the correction equation
      write(*,'(a,a)') 'Calling b2uxus from ', name
      call b2uxus (nCv, mpg, aa, itcnt, resmb, corub, name)           !srv 18.03.10

      if (iout.ne.0) then                                             !srv 17.06.02 {
#ifdef WG_TODO
       call my_out_us(70,nCv,0,aa(-1,-1,-4),name//'_aam4')
       call my_out_us(70,nCv,0,aa(-1,-1,-3),name//'_aam3')
       call my_out_us(70,nCv,0,aa(-1,-1,-2),name//'_aam2')
       call my_out_us(70,nCv,0,aa(-1,-1,-1),name//'_aam1')
       call my_out_us(70,nCv,0,aa(-1,-1, 0),name//'_aam0')
       call my_out_us(70,nCv,0,aa(-1,-1, 1),name//'_aap1')
       call my_out_us(70,nCv,0,aa(-1,-1, 2),name//'_aap2')
       call my_out_us(70,nCv,0,aa(-1,-1, 3),name//'_aap3')
       call my_out_us(70,nCv,0,aa(-1,-1, 4),name//'_aap4')
#endif
       call my_out_us(70,nCv,0,resmb,name//'_resmb')
       call my_out_us(70,nCv,0,corub,name//'_corub')
       call my_out_us(70,nCv,0,pccb(1,0),name//'_pccbx')
       call my_out_us(70,nCv,0,pccb(1,1),name//'_pccby')
      endif                                                           !srv 17.06.02 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2usmo

      end
