*-----------------------------------------------------------------------
*.specification

      subroutine b2tqin (nCv, ns, nscx, iscx, switch, geo,
     &                   pl, rt, sigin)
      use b2mod_types
      use b2mod_math
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_plasma
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2tqin
#ifdef B25_EIRENE
      use eirmod_wneutrals
     & , only : dab2
#else
      use b2mod_eirdiag
     & , only : dab2
#endif
      use b2mod_neutrals_namelist
     &  , only : b2eatcr
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, ns, nscx
      integer iscx(0:nscx-1)
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (B2Plasma), intent (in) :: pl
      type (B2Rates), intent (in) :: rt
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  sigin(nCv,0:1,0:ns-1,0:nscx-1)                                 !sv 01.07.99 !xpb
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TQIN computes the classical ion-neutral perpendicular conductivity.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the ion-neutral perpendicular conductivity
*     sigin: electrical conductivity.
*     This coefficient occur in the equation as indicated below.
*     Simplifed electric current.
*       fchx = sigin*ehx
*       fchy = sigin*ehy
*     Here, ehx and ehy are the components of the ... modified electric
*     field: ehx = -(d/dx)po - (1/(qe*ni))*(d/dx)(ni*ti)+B*hx*Vyn,
*     ehy = -(d/dy)po - (1/(ni*te))*(d/dy)(ne*te)-B*hy*Vpn.
*     See ...
*
*     The equations indicated above do not describe all the relevant
*     transport properties of a plasma. The most significant omission
*     is the absence of diamagnetic transport terms.
*
*
*  4. references
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, is, is0, k
      real (kind=R8) ::
     *   rcx, t0, tav
      real (kind=R8) :: na0_kin(nCv), na0_fluid(nCv), na0_tot(nCv) !! fluid and kinetic densities of main hydrogen neutral atom
*   ..procedures
      intrinsic sqrt, log
      external xertst, sfill
      external b2xvsg
*   ..initialisation

*-----------------------------------------------------------------------
*.documentation-internal
*
*     Description of some local variables:
*
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tqin')
*   ..test nCv, ns
      call xertst (0.lt.nCv, 'faulty argument nCv')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
      if (ncall_b2tqin.lt.3) then
*    ..test bb
       call b2xvsg (nCv, geo%cvBb(1,3), 1, 'bb3', '.gt.')
*    ..test state
       call b2xvsg (nCv, pl%ti, 1, 'ti', '.gt.')
       call b2xvsg (nCv*ns, pl%na, 1, 'na', '.gt.')
      endif

* ..compute classical sigin
*   ..compute sigin.
      call sfill (nCv*2*ns*nscx, 0.0_R8, sigin, 1)
*   ..we use the charge exchange rate as computed from b2frates
*   ..this is only nonzero for hydrogen neutral on hydrogen ion CX
#if (defined(_OPENMP) && !defined(NO_OPENMP_B2TQIN))
!svdk parallellize over charge exchange species 
!$omp parallel do default(none)
!$omp& private(k, is0, is, iCv, tav, t0, rcx, na0_kin, na0_fluid)
!$omp& private(na0_tot)
!$omp& shared(nCv, is_neutral, zn, ns, nscx, iscx, dab2, b2eatcr)
!$omp& shared(pl , am, switch, sigin, geo, rt)
#endif
      do k = 0, nscx-1
        is0 = iscx(k)
        if (switch%use_eirene.ne.0) then
          na0_kin = dab2(:,b2eatcr(is0),1)
        else
          na0_kin = 0.0_R8
        endif
        na0_fluid = pl%na(:,is0)
        na0_tot = na0_kin + na0_fluid !! sum for hybrid fluid-kinetic simulations

        do is = 0, ns-1
         if (.not.is_neutral(is).and.nint(zn(is)).eq.1) then
          do iCv = 1, nCv
            tav = (pl%ti(iCv)+pl%tn(iCv))/2.0_R8
            t0 = log(tav/(am(is0)*ev))
            rcx = switch%b2sqcx_phm0*
     &       expu(rt%rlcx(iCv,0,is,k)+rt%rlcx(iCv,1,is,k)*t0)
            if (switch%b2tqin_csigin_style.eq.0) then
              ! From 5.0/5.1: Using the CX rate
              if (switch%transport_afn.eq.0) then !! default
                sigin(iCv,0,is,k) = (pl%na(iCv,is)*am(is)*mp)
     &          *rcx*na0_tot(iCv)/(2.0_R8*geo%cvBb(iCv,3)**2)
              else
!! in case of advanced fluid neutrals, rcx = Kcxm and in that case the factor 1/2 should not be there
!! ( see "Van Uytven, Wim, et al. "Effect of drifts and currents on the validity of a fluid model for the atoms in the plasma edge." Nuclear Materials and Energy 33 (2022): 101255.")
                sigin(iCv,0,is,k) = (pl%na(iCv,is)*am(is)*mp)
     &          *rcx*na0_tot(iCv)/(geo%cvBb(iCv,3)**2)
              endif
            else
              ! From 5.2: Using the diffusivity of room temperature hydrogen gas
              sigin(iCv,0,is,k) = (pl%na(iCv,is)*am(is)*mp)      !srv 01.07.99
     &         *(pl%na(iCv,is0)*3.2e-15_R8)
     &         *sqrt(pl%ti(iCv)/ev/0.026_R8)/(2.0_R8*geo%cvBb(iCv,3)**2)
            endif
            sigin(iCv,1,is,k) = sigin(iCv,0,is,k)
          enddo
         endif
        enddo
      enddo
#if (defined(_OPENMP) && !defined(NO_OPENMP_B2TQIN))
!$omp end parallel do
#endif

* ..return
      ncall_b2tqin = ncall_b2tqin + 1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tqin

      end subroutine b2tqin
