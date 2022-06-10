*-----------------------------------------------------------------------
*.specification

      subroutine b2tqca (nCv, ns, switch, geo, pl, dv, rt, st_ext, 
     &                   vsa, hci)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_plasma
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2tqca, ncall_b2tlnl
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, ns
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (in) :: dv
      type (B2Rates), intent (in) :: rt
      type (B2StateExt), intent (in) :: st_ext
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  vsa(nCv,0:1,0:ns-1), hci(nCv,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TQCA computes the classical contributions to the atom
*     transport coefficients.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the classical contributions to the following
*     atom transport coefficients.
*     - vsa: viscosity.
*     - hci: atom thermal conductivity.
*     These coefficients occur in the equations as indicated below.
*     Parallel momentum transport.
*       momx(is) = -vsax(is)*(d/dx)up(is)
*       momy(is) = -vsay(is)*(d/dy)up(is)
*     Atom heat flux.
*       fhix = -hcix*(d/dx)ti
*       fhiy = -hciy*(d/dy)ti
*     The particle and parallel momentum transport is a property of
*     each atomic species, the other transport laws involve summation
*     over all atomic species.
*
*     The equations indicated above do not describe all the relevant
*     transport properties of a plasma. The most significant omission
*     is the absence of diamagnetic transport terms. Also, one might
*     define temperatures and heat fluxes for each individual atom
*     species. These tasks are left for future work.
*
*
*  4. references
*
*     /1/ R. Balescu: Transport Processes in Plasmas, Vol. 1.
*     North Holland, Amsterdam, 1988, pp. 230-239.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, is, js
      character chns*3                                                !srv 16.06.08
      real (kind=R8) :: 
     *  fki, fxi, eta, ctaup, ctaui, ci, 
     *  t, t0, t1, tau(nCv), tnp, w                             !srv 11.09.09
      parameter (fki=2.253_R8, eta=0.96_R8, ci = 3.9e0_R8)
*   ..procedures
      intrinsic abs, sqrt, log, max
      real (kind=R8) ::
     *  pitx
      external subini, subend, xertst, sfill
      external b2xvsg
      pitx(iCv) = abs(geo%cvBb(iCv,0))/geo%cvBb(iCv,3)
*   ..initialisation

*-----------------------------------------------------------------------
*.documentation-internal
*
*     Description of some local variables:
*
*     lnlam is the log of the plasma parameter.
*     fki is the dimensionless parallel atom thermal conductivity.
*     eta is the dimensionless parallel atom viscosity.
*     (see table 4.1 of Balescu for fki and eta.)
*
*     Description of some local procedures:
*
*     pitx(iCv) is the local pitch of the magnetic field.
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tqca')
*   ..test nCv, ns
      call xertst (0.le.nCv, 'faulty argument nCv')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
      if (ncall_b2tqca.lt.3) then
*    ..test bb
        call b2xvsg (nCv, geo%cvBb(1,3), 1, 'bb3', '.gt.')
*    ..test state
        call b2xvsg (nCv*ns, pl%na, 1, 'na', '.gt.')
        call b2xvsg (nCv, pl%ti, 1, 'ti', '.gt.')
        call b2xvsg (nCv, dv%ne2, 1, 'ne2', '.gt.')
      endif
*   ..compute ctaup
*     Following Balescu; Braginskii has ctaup different by a factor
*     sqrt(2). The mathematical expression is:
*       t0 = (3/(4*sqrt(2*pi)))*(4*pi*eps0/qe**2)**2*ev*sqrt(ev)
*       ctaup = t0*sqrt(mp)/lnlam
*     The numerical value in SI units is: ctaup = 1.4743e13/lnlam.
*     The following calculation of t0 is designed to avoid overflow
*     even for IEEE single precision arithmetic.
!xpb The lnlam term has been removed and is put back in later so it may
!xpb vary with space and species or be set by user.
      t1 = (4.0_R8*pi*eps0/qe)/qe
      t0 = (3.0_R8/(4.0_R8*sqrt(2.0_R8*pi)))*t1*(t1*ev*sqrt(ev))
      ctaup = t0*sqrt(mp)						 !xpb

*   ..compute the Coulomb logarithm
      call b2tlnl (nCv, switch, switch%icase_ii, pl%te, pl%ti, 
     &                                              dv%ne, dv%lnlam)      !srv 01.07.09  20.09.11
* ..compute vsa, hci
*   ..compute vsax, hcix
*     (initialise hcix to 0)
      call sfill (nCv, 0.0_R8, hci(1,0), 1)
      do is = 0, ns-1
       if (is_neutral(is)) then
*       (neutral vsax and hcix are strictly anomalous)
        call sfill (nCv, 0.0_R8, vsa(1,0,is), 1)
       else
        do iCv = 1, nCv
          t0 = pl%ti(iCv)/ev
          tau(iCv) = (ctaup/dv%lnlam(iCv))*sqrt(am(is))*t0*sqrt(t0)/          !xpb !srv 01.07.09
     &      (rt%rz2(iCv,is)*dv%ne2(iCv))
          tnp = tau(iCv)*pl%na(iCv,is)*pitx(iCv)**2
          if (switch%b2tqca_model.eq.1) then  ! SOLPS5.0 model : Balescu 29-moment result
            hci(iCv,0) = hci(iCv,0)+
     &       2.5_R8*(pl%ti(iCv)/(mp*am(is)))*fki*tnp
c DPC/XPB changed to 4/3*  at the strong recommendation of Rozhansky
c WDK/SRV 21.07.2016 factor 4/3 removed from vsa, but added explicitly to cvsahz (not cvsa),
c         according to old model.eq.3 treatment
            vsa(iCv,0,is) = 
     &       switch%b2tqca_phm0*
     &       eta*sqrt(2.0_R8)*pl%ti(iCv)*tnp
          else if (switch%b2tqca_model.eq.2) then ! old SOLPS4.0 model : Braginskii
            w = 0.0_R8
            t = 0.0_R8
            do js = 0, ns-1
              w = w + rt%rz2(iCv,js) * pl%na(iCv,js) * 
     &            sqrt(2.0_R8*mp*am(is)*am(js)/(am(is)+am(js))) 
              t = t + rt%rz2(iCv,js) * pl%na(iCv,js) / 
     &            sqrt((am(is)+am(js))/2.0_R8) 
            enddo
            do js = 0, st_ext%ns-1
              w = w + st_ext%za2(iCv,js) * st_ext%na(iCv,js) *
     &            sqrt(2.0_R8*mp*am(is)*st_ext%am(js)/
     &            (am(is)+st_ext%am(js)))
              t = t + st_ext%za2(iCv,js) * st_ext%na(iCv,js) /
     &            sqrt((am(is)+st_ext%am(js))/2.0_R8)
            enddo
cxpb            ctaui = ctaup*sqrt(2.0_R8)/lnlam/rz2(iCv,is)  ! lower precision as per 4.0
	    ctaui = 2.1e13_R8 / (dv%lnlam(iCv)*rt%rz2(iCv,is))             !srv 01.07.09
	    fxi = ci * ctaui / (ev*sqrt(ev*mp))
            hci(iCv,0) = hci(iCv,0) + fxi*(pl%na(iCv,is)/w)*
     &       pitx(iCv)**2*(pl%ti(iCv)**2*sqrt(pl%ti(iCv)))
cwdk      no factor 4/3 here, for consistent treatment with model.eq.1
            vsa(iCv,0,is) = 
     &       switch%b2tqca_phm0*
     &       eta*ctaui/(ev*sqrt(ev))*(pl%na(iCv,is)/t)*
     &       pitx(iCv)**2*(pl%ti(iCv)**2*sqrt(pl%ti(iCv)))
          else
            call xerrab ('b2tqca--invalid value b2tqca_model > 2. ')
          endif
        enddo
        if (switch%b2npht_iout.ne.0
     &     .or.switch%b2tqca_iout.ne.0.or.switch%iout_b2wdat.eq.4) then
         write (chns,'(i3.3)') is
         call my_out_us(70,nCv,0,tau,'b2tqca_taua'//chns)
	endif
       endif
      enddo
*   ..compute vsay, hciy
*!!   (classical perpendicular transport is set to 0)
      do is = 0, ns-1
       call sfill (nCv, 0.0_R8, vsa(1,1,is), 1)
      enddo
      call sfill (nCv, 0.0_R8, hci(1,1), 1)

      if(switch%b2tqca_iout.ne.0 .or. switch%b2npmo_iout.ne.0) then       !srv 16.06.08 30.01.10 {
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        call my_out_us(70,nCv,0,vsa(1,0,is),'b2tqca_vsax'//chns)
        call my_out_us(70,nCv,0,vsa(1,1,is),'b2tqca_vsay'//chns)
!        call my_out_us(70,nCv,0,pl%na(1,is),'b2tqca_na'//chns)
       enddo
!        call my_out_us(70,nCv,0,lnlam,'b2tqca_lnlam')
!        call my_out_us(70,nCv,0,pl%ti,'b2tqca_ti')
!        call my_out_us(70,nCv,0,dv%ne,'b2tqca_ne')
!        call my_out_us(70,nCv,0,dv%ne2,'b2tqca_ne2')
!        call my_out_us(70,nCv,0,geo%cvBb(1,0),'b2tqca_Bx')
!        call my_out_us(70,nCv,0,geo%cvBb(1,3),'b2tqca_BB')
      endif                                                          !srv 16.06.08 }

      if(switch%b2tqca_iout.ne.0 .or. switch%iout_b2wdat.eq.4) then  !srv 16.06.08 30.01.10 {
        call my_out_us(70,nCv,0,dv%lnlam,'b2tqca_lnlam')
      endif

      if(switch%b2tqca_iout.ne.0 .or. switch%b2npht_iout.ne.0) then
        call my_out_us(70,nCv,0,hci(1,0),'b2tqca_hcix')
        call my_out_us(70,nCv,0,hci(1,1),'b2tqca_hciy')
      endif                                                          !srv 16.06.08 30.01.10 }

* ..return
      ncall_b2tqca = ncall_b2tqca+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tqca

      end subroutine b2tqca
