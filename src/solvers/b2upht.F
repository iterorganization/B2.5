*-----------------------------------------------------------------------
*.specification

      subroutine b2upht (nCv, nFc, nVx, ns, switch, geo, mpg, nregionv, 
     &  po_solve, ua_solve, solvereg, solvireg, solvnreg, solvtreg,       !srv 22.05.18
     &  solvpreg, solvmreg, rxf, pcm0, cortt, corte, corti, cortn,
     &  pccm0, na, ua, te, ti, tn, po, ne, ni)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, nregionv
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *  rxf, pcm0, cortt(nCv), corte(nCv),
     &  corti(nCv), cortn(nCv), pccm0(nCv)
      logical po_solve, ua_solve,
     &  solvereg(0:nregionv), solvireg(0:nregionv),
     &  solvnreg(0:nregionv), solvtreg(0:nregionv),                             !srv 22.05.18
     &  solvpreg(0:nregionv), solvmreg(0:ns-1,0:nregionv)
*   ..input/output arguments
      real (kind=R8) ::
     *  na(nCv,0:ns-1), ua(nCv,0:ns-1),
     &  te(nCv), ti(nCv), tn(nCv), po(nCv),
     &  ne(nCv), ni(nCv,0:1)
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
      integer ncall, iCv, is, ireg
      real (kind=R8) ::
     *  ctt, dte(nCv), dti(nCv), dtn(nCv), t0, t1, t2,
     *  wrk(nCv), wrk0(nCv), wrkv(nVx)
      save ncall
*   ..procedures
      intrinsic min, max
      real (kind=R8) ::
     *  med
      external subini, subend, xertst, b2xvsg, xerrab
      med(t0,t1,t2) = min(max(t0,t1),max(t1,t2),max(t2,t0))
*   ..initialisation
      data ncall/0/


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
      if(ncall.lt.3) then
        call b2xvsg (nCv*ns, na, 1, 'na', '.gt.')
        call b2xvsg (nCv*2, ni, 1, 'ni', '.gt.')
        call b2xvsg (nCv, geo%cvVol, 1, 'vol', '.gt.')
        call b2xvsg (nCv, te, 1, 'te', '.gt.')
        call b2xvsg (nCv, ti, 1, 'ti', '.gt.')
        call b2xvsg (nCv, tn, 1, 'tn', '.gt.')
        call b2xvsg (nCv, ne, 1, 'ne', '.gt.')
      endif 

* ..apply the correction
*   ..adjust te, ti, po
      do iCv = 1, nCv
        ireg = mpg%cvReg(iCv)
*!!   (stylec is a numerical switch, needs experiments)
        if (switch%stylec.eq.0) then
         if (solvereg(ireg).and.solvireg(ireg).and.solvtreg(ireg)) then     !srv 22.05.18
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

         if (solvereg(ireg).and.solvireg(ireg).and.solvtreg(ireg)) then     !srv 22.05.18
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
        if (solvereg(ireg)) te(iCv) = max(1.0e-5_R8*ev,te(iCv)+dte(iCv))	 !xpb
        if ((solvpreg(ireg).and.po_solve).and.switch%pot_eq.eq.1.and.
     &      iCv.le.mpg%nCi) then
          po(iCv) = po(iCv)+dte(iCv)/qe
        else if (switch%pot_eq.eq.2) then
          po(iCv) = 3.1_R8*te(iCv)/qe
        end if
        if (solvireg(ireg)) ti(iCv) = max(1.0e-5_R8*ev,ti(iCv)+dti(iCv))	 !xpb
        if (solvnreg(ireg)) tn(iCv) = max(1.0e-5_R8*ev,tn(iCv)+dtn(iCv))	 
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

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2upht

      end subroutine b2upht
