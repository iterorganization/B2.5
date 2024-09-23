*-----------------------------------------------------------------------
*.specification

      subroutine b2usmo (nCv, nFc, nVx, switch, geo, mpg,
     &  nregionv, solvereg, itcnt, rxg, rob, rzb, pb, pz, ub,        !srv 09.01.01
     &  smb, flcb, cvsb, resmb, ctcfb, corub, pccb, aa, aad0, name)  !srv 06.04.07
      use b2mod_types
      use b2mod_switches
      use b2us_geo
      use b2us_map
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2usmo
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, nregionv, itcnt
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      character*(*) name                                             !srv 18.05.02
      real (kind=R8) ::
     *  rxg, rob(nCv), rzb(nCv), pb(nCv), pz(nCv),
     &  ub(nCv), smb(nCv,0:3), flcb(nFc,0:1), cvsb(nFc,0:1),
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
*     B2USMO computes a correction towards the solution of the
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
      integer iCv, iFc, iCv1, iCv2, iVx1, iVx2, inCv, ireg
      real (kind=R8) ::
     *  t1, t2, t3, csb, csb1,
     *  t11, t22
      logical lnonzero
*   ..procedures
      external xertst, ipgetr, sfill
      external b2xvsg, b2ux5p, intfaceh, intfacev
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2usmo')
*   ..test nCv, nFc
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
*   ..test itcnt
      call xertst (0.le.itcnt, 'faulty argument itcnt')
*   ..extensive tests on first few calls
      if (ncall_b2usmo.lt.3) then
*    ..test sign of rob, rzb, pb, pz, smb1, smb3, cvsb
       call b2xvsg (nCv, rob, 1, 'rob', '.gt.')
       call b2xvsg (nCv, rzb, 1, 'rzb', '.gt.')
       call b2xvsg (nCv, pb, 1, 'pb', '.gt.')
       call b2xvsg (nCv, pz, 1, 'pz', '.gt.')
       call b2xvsg (nCv, smb(1,1), 1, 'smb1', '.le.')
       call b2xvsg (nCv, smb(1,3), 1, 'smb3', '.le.')
      endif

* ..apply discretization scheme
      call calccoef(nCv,nFc,nVx,switch%b2npmo_discr_meth,
     &                                geo,flcb,cvsb,flcb0,cvsb0)

* ..compute the correction
      aa = 0.0_R8
*   ..compute the nine-point matrix
! $AD II-LOOP
      do iCv = 1, nCv
        ireg = mpg%cvReg(iCv)
        if (solvereg(ireg)) then
         csb = sqrt(pz(iCv)/rzb(iCv))
         csb1 = sqrt(ub(iCv)**2+csb**2)
         aad0(iCv) = -smb(iCv,1)-
     &    smb(iCv,3)*rob(iCv)+
     &    switch%cfc0*geo%cvVol(iCv)*geo%cvHz(iCv)/geo%cvHx(iCv)*       !srv 09.01.01
     &    geo%cvBb(iCv,0)/geo%cvBb(iCv,3)*csb*rob(iCv)
         aa(mpg%cvNvP(iCv,1)) = aad0(iCv)
         if (switch%b2mndt_style.lt.2) then
           aa(mpg%cvNvP(iCv,1)) = aad0(iCv) +abs(resmb(iCv))/(rxg*csb1)
         else
           aa(mpg%cvNvP(iCv,1)) = aad0(iCv)
         endif
        else
         aa(mpg%cvNvP(iCv,1)) = 1.0_R8
         aad0(iCv) = 0.0_R8
        endif
      enddo
! $AD II-LOOP
      do iFc = 1, nFc
        iCv1 = mpg%fcCv(iFc,1) ! cell number of c1
        iCv2 = mpg%fcCv(iFc,2) ! cell number of c2
*     ..contributions from flcb0(:,0), cvsb0(:,0) ('5-point' contributions)
        t1 = - 0.5_R8*flcb0(iFc,0) - cvsb0(iFc,0)
        t2 =   0.5_R8*flcb0(iFc,0) - cvsb0(iFc,0)
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
*     ..contributions from flcb0(:,1), cvsb0(:,1) ('9-point' contributions)
        iVx1 = mpg%fcVx(iFc,1) ! number of v1
        iVx2 = mpg%fcVx(iFc,2) ! number of v2
        t1 = - 0.5_R8*flcb0(iFc,1) - cvsb0(iFc,1)
        t2 =   0.5_R8*flcb0(iFc,1) - cvsb0(iFc,1)
        t11 = 0.0_R8
        do inCv = 1, mpg%vxCvP(iVx1,2)
          t11 = t11 + 1.0_R8/geo%vxVol(mpg%vxCvP(iVx1,1)+inCv-1)
        enddo
        t22 = 0.0_R8
        do inCv = 1, mpg%vxCvP(iVx2,2)
          t22 = t22 + 1.0_R8/geo%vxVol(mpg%vxCvP(iVx2,1)+inCv-1)
        enddo
        ! contribution to equation of iCv1
        if (solvereg(mpg%cvReg(iCv1))) then
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
        endif
        ! contribution to equation of iCv2
        if (solvereg(mpg%cvReg(iCv2))) then
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
        endif
      enddo
*   ..compute pccb
      pccb = 0.0_R8
! $AD II-LOOP
      do iCv = 1, nCv
        ireg = mpg%cvReg(iCv)
        if (solvereg(ireg)) then
         pccb(iCv,0) = geo%cvVol(iCv)*geo%cvHz(iCv)*                    !srv 09.01.01
     .    geo%cvBb(iCv,0)/(geo%cvBb(iCv,3)*aa(mpg%cvNvP(iCv,1)))   !srv 11.01.13
         pccb(iCv,1) = geo%cvVol(iCv)*geo%cvHz(iCv)*                    !srv 09.01.01
     .    geo%cvBb(iCv,0)/(geo%cvBb(iCv,3)*aad0(iCv))              !srv 11.01.13
        else
         resmb(iCv) = 0.0_R8   ! take opportunity to correct RHS
        endif
      enddo
*   ..contribution from linearization of centrifugal force
! $AD II-LOOP
      do iCv = 1, nCv
        aa(mpg%cvNvP(iCv,1)) = aa(mpg%cvNvP(iCv,1)) -
     &    min(0.0_R8,ctcfb(iCv))
      enddo
! $AD II-LOOP
      do iCv = 1, nCv
        lnonzero=.false.
        do inCv = 1,mpg%cvNvP(iCv,2)
          lnonzero = lnonzero.or.(aa(mpg%cvNvP(iCv,1)+inCv-1).ne.0.0_R8)
        enddo
        if (.not.lnonzero) then
          aa(mpg%cvNvP(iCv,1)) = 1.0_R8
          resmb(iCv) = 0.0_R8
        endif
      enddo

*   ..solve the correction equation
      write(*,'(a,a)') 'Calling b2uxus from ', name
      call b2uxus (nCv, mpg, aa, itcnt, resmb, corub, name)           !srv 18.03.10

      if (switch%b2usmo_iout.ne.0) then                               !srv 17.06.02 {
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
       call my_out_us(70,nCv,0,rob,name//'_rob')
       call my_out_us(70,nCv,0,pb,name//'_pb')
       call my_out_us(70,nCv,0,smb(:,0),name//'_smb0')
       call my_out_us(70,nCv,0,smb(:,1),name//'_smb1')
       call my_out_us(70,nCv,0,smb(:,2),name//'_smb2')
       call my_out_us(70,nCv,0,smb(:,3),name//'_smb3')
       call my_out_us(70,nFc,1,flcb,name//'_flcb')
       call my_out_us(70,nFc,1,cvsb,name//'_cvsb')

      endif                                                           !srv 17.06.02 }

* ..return
      ncall_b2usmo = ncall_b2usmo+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2usmo

      end subroutine b2usmo
