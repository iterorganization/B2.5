*-----------------------------------------------------------------------
*.specification

      subroutine b2upht (nCv, nFc, nVx, ns, switch, geo, mpg, nregionv, 
     &  po_solve, ua_solve, solvereg, solvireg, solvnreg, solvtreg,       !srv 22.05.18
     &  solvpreg, solvmreg, solvkreg, solvzreg, rxf, pcm0, cortt, corte,
     &  corti, cortn, corkt, corzt, pccm0, na, ua, te, ti, tn, kt, zt,
     &  po, ne, ni)
      use b2mod_types
      use b2mod_math
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_numerics_namelist
     , , only : corr_core_dt, add_te_corr_to_po              ! IYS 19.12.2017, 21.09.2018
      use b2mod_switches
      use b2us_geo
      use b2us_map
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2upht
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, nregionv
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *  rxf, pcm0, cortt(nCv), corte(nCv),
     &  corti(nCv), cortn(nCv), corkt(nCv), corzt(nCv), pccm0(nCv)
      logical po_solve, ua_solve,
     &  solvereg(0:nregionv), solvireg(0:nregionv),
     &  solvnreg(0:nregionv), solvtreg(0:nregionv),                             !srv 22.05.18
     &  solvpreg(0:nregionv), solvmreg(0:ns-1,0:nregionv),
     &  solvkreg(0:nregionv), solvzreg(0:nregionv)
*   ..input/output arguments
      real (kind=R8) ::
     *  na(nCv,0:ns-1), ua(nCv,0:ns-1),
     &  te(nCv), ti(nCv), tn(nCv), po(nCv),
     &  ne(nCv), ni(nCv,0:1), kt(nCv), zt(nCv)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2UPHT adjusts the state of the plasma towards the solution
*     of the heat equations.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, is, ireg, iFt
      real (kind=R8) ::
     *  ctt, dte(nCv), dti(nCv), dtn(nCv), dkt(nCv), dzt(nCv), t0,
     *  wrk(nCv), wrk0(nCv), wrkv(nVx)
      real (kind=R8) ::
     *  dte1(nCv), dti1(nCv),                                       !lk 26.10.17
     *  tt0(mpg%nFt), tt1(mpg%nFt), tt2(mpg%nFt), tt3(mpg%nFt)      !lk 26.10.17
*   ..procedures
      intrinsic min, max
      external xertst, b2xvsg, xerrab
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2upht')
*   ..test nCv, nFc, ns
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test rxf, pcm0
      call xertst (0.0_R8.le.rxf.and.rxf.le.1.0_R8,
     & 'faulty argument rxf')
      call xertst (0.0_R8.le.pcm0, 'faulty argument pcm0')
*   ..test sign of na, ni, vol, te, ti, ne
      if(ncall_b2upht.lt.3) then
        call b2xvsg (nCv*ns, na, 1, 'na', '.gt.')
        call b2xvsg (nCv*2, ni, 1, 'ni', '.gt.')
        call b2xvsg (nCv, geo%cvVol, 1, 'vol', '.gt.')
        call b2xvsg (nCv, te, 1, 'te', '.gt.')
        call b2xvsg (nCv, ti, 1, 'ti', '.gt.')
        call b2xvsg (nCv, tn, 1, 'tn', '.gt.')
        call b2xvsg (nCv, ne, 1, 'ne', '.gt.')
        call b2xvsg (nCv, kt, 1, 'kt', '.ge.')
        call b2xvsg (nCv, zt, 1, 'zt', '.ge.')
      endif

* ..apply the correction
*   ..adjust te, ti, po
      tt0 = 0.0e0_R8
      tt1 = 0.0e0_R8
      tt2 = 0.0e0_R8
      tt3 = 0.0e0_R8
      do iCv = 1, nCv
        ireg = mpg%cvReg(iCv)
*!!   (stylec is a numerical switch, needs experiments)
        if (switch%stylec.eq.0) then
         if (solvereg(ireg).and.solvireg(ireg).and.solvtreg(ireg)) then  !srv 22.05.18
          dte(iCv) = rxf*max(-te(iCv),(corte(iCv)+cortt(iCv))/2.0_R8)
          dti(iCv) = rxf*max(-ti(iCv),(corti(iCv)+cortt(iCv))/2.0_R8)
         else
          dte(iCv) = rxf*max(-te(iCv),corte(iCv))
          dti(iCv) = rxf*max(-ti(iCv),corti(iCv))
         endif
         if (solvereg(ireg).and.solvireg(ireg).and.solvnreg(ireg)) then
          dtn(iCv) = rxf*max(-tn(iCv),(cortn(iCv)+cortt(iCv))/2.0_R8)
         elseif (solvnreg(ireg)) then
          dtn(iCv) = rxf*max(-tn(iCv),cortn(iCv))
         else 
          dtn(iCv) = 0.0_R8
         endif
        else if (switch%stylec.eq.1) then

         if (solvereg(ireg).and.solvireg(ireg).and.solvtreg(ireg)) then  !srv 22.05.18
          ctt = cortt(iCv) - med(corte(iCv),cortt(iCv),corti(iCv))
          dte(iCv) = rxf*max(-te(iCv),corte(iCv)+ctt)
          dti(iCv) = rxf*max(-ti(iCv),corti(iCv)+ctt)
         else
          dte(iCv) = rxf*max(-te(iCv),corte(iCv))
          dti(iCv) = rxf*max(-ti(iCv),corti(iCv))
         endif
        else if (switch%stylec.eq.2) then
*        !! This for testing purposes only
         dte(iCv) = rxf*max(-te(iCv),corte(iCv))
         dti(iCv) = 0.0_R8
        else if (switch%stylec.eq.3) then
*        !! This for testing purposes only
         dte(iCv) = 0.0_R8
         dti(iCv) = rxf*max(-ti(iCv),corti(iCv))
        else if (switch%stylec.eq.4) then
*        !! This for testing purposes only
         dte(iCv) = rxf*max(-te(iCv),corte(iCv))
         dti(iCv) = rxf*max(-ti(iCv),corti(iCv))
        endif
        dkt(iCv)  = rxf*max(-kt(iCv),corkt(iCv))
        dzt(iCv)  = rxf*max(-zt(iCv),corzt(iCv))
        dte1(iCv) = dte(iCv)
        dti1(iCv) = dti(iCv)
*   ..smooth the correction for inner flux surfaces                 !lk 26.10.2017{
        if(mpg%cvOnClosedSurface(iCv).and.(iCv.le.mpg%nCi)) then    ! IYS 27.03.2019
          iFt = mpg%cvFt(iCv)
          tt0(iFt) = tt0(iFt) + geo%cvVol(iCv)*ni(iCv,1)
          tt1(iFt) = tt1(iFt) + geo%cvVol(iCv)*ni(iCv,1)*dti1(iCv)
          tt2(iFt) = tt2(iFt) + geo%cvVol(iCv)*ne(iCv)
          tt3(iFt) = tt3(iFt) + geo%cvVol(iCv)*ne(iCv)*dte1(iCv)
        endif
      enddo

      do iCv = 1, nCv
        ireg = mpg%cvReg(iCv)
        if(mpg%cvOnClosedSurface(iCv).and.(iCv.le.mpg%nCi)) then    ! IYS 27.03.2019
          iFt = mpg%cvFt(iCv)
          dti1(iCv) = tt1(iFt)/tt0(iFt) +
     *         corr_core_dt*(dti1(iCv)-tt1(iFt)/tt0(iFt))           ! IYS 19.12.2017
          dte1(iCv) = tt3(iFt)/tt2(iFt) +
     *         corr_core_dt*(dte1(iCv)-tt3(iFt)/tt2(iFt))           ! IYS 19.12.2017
        endif
        dte(iCv) = dte1(iCv)
        dti(iCv) = dti1(iCv)                                        !srv 25.07.17
        if (solvereg(ireg)) te(iCv) = max(1.0e-5_R8*ev,te(iCv)+dte(iCv))  !xpb
        if ((solvpreg(ireg).and.po_solve).and.switch%pot_eq.eq.1.and.
     &      iCv.le.mpg%nCi.and.add_te_corr_to_po(ireg)) then
          po(iCv) = po(iCv)+dte(iCv)/qe
        else if (switch%pot_eq.eq.2) then
          po(iCv) = 3.1_R8*te(iCv)/qe
        end if
        if (solvireg(ireg)) ti(iCv) = max(1.0e-5_R8*ev,ti(iCv)+dti(iCv))  !xpb
        if (solvnreg(ireg)) tn(iCv) = max(1.0e-5_R8*ev,tn(iCv)+dtn(iCv))
        if (switch%solve_keps.gt.0.and.solvkreg(ireg))
     &    kt(iCv) = max(1.0e-5_R8*ev,kt(iCv)+dkt(iCv))
        if (switch%solve_keps.gt.1.and.solvzreg(ireg))
     &    zt(iCv) = max(1.0e-8_R8*ev,zt(iCv)+dzt(iCv))
      enddo
*   ..adjust ua
      wrk0 = ne*dte + ni(:,0)*dti
      call gradc_p(nCv, nFc, nVx, 0, geo, mpg, wrk0, wrkv, wrk)
      if (ua_solve) then
       do iCv = 1, mpg%nCi
        ireg = mpg%cvReg(iCv)
        t0 = pcm0*pccm0(iCv)*wrk(iCv)
        do is = 0, ns-1
         if (solvmreg(is,ireg)) ua(iCv,is) = ua(iCv,is)-t0
        enddo
       enddo
      endif

      if (switch%tn_style.le.1) then
        tn = ti
      endif

* ..return
      ncall_b2upht = ncall_b2upht+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2upht

      end subroutine b2upht
