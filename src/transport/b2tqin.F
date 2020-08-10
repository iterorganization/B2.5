*-----------------------------------------------------------------------
*.specification

      subroutine b2tqin (nCv, ns, nscx, iscx, geo, pl, rt, sigin)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_plasma
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, ns, nscx
      integer iscx(0:nscx-1)
      type (geometry), intent (in) :: geo
      type (B2Plasma), intent (in) :: pl
      type (B2Rates), intent (in) :: rt
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *	sigin(nCv,0:1,0:ns-1,0:nscx-1)                                      !sv 01.07.99 !xpb
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
      integer ncall, iCv, is, is0, k
      real (kind=R8) ::
     *   phm0, rcx, expu, t0, tav
      save ncall, phm0
*   ..procedures
      intrinsic sqrt, log, exp, max, dim
      external subini, subend, xertst, sfill, ipgetr, expu
      external b2xvsg
*   ..initialisation
      data ncall/0/, phm0/1.0_R8/

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
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2sqcx_phm0', phm0)
       call xertst (0.0_R8.le.phm0, 'faulty internal parameter phm0')
      endif
*   ..test nCv, ns
      call xertst (0.le.nCv, 'faulty argument nCv')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
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
      do k = 0, nscx-1
        is0 = iscx(k)
        do is = 0, ns-1
         if (.not.is_neutral(is).and.nint(zn(is)).eq.1) then
          do iCv = 1, nCv
            tav = (pl%ti(iCv)+pl%tn(iCv))/2.0_R8
            t0 = log(tav/(am(is0)*ev))
            rcx = phm0*
     &       expu(rt%rlcx(iCv,0,is,k)+rt%rlcx(iCv,1,is,k)*t0)
! From 5.2: Using the diffusivity of room temperature hydrogen gas
            sigin(iCv,0,is,k) = (pl%na(iCv,is)*am(is)*mp)      !srv 01.07.99
     &       *(pl%na(iCv,is0)*3.2e-15_R8)
     &       *sqrt(pl%ti(iCv)/ev/0.026_R8)/(2.0_R8*geo%cvBb(iCv,3)**2)
! From 5.0/5.1: Using the CX rate
!            sigin(iCv,0,is,k) = (pl%na(iCv,is)*am(is)*mp)
!     &       *rcx*pl%na(iCv,is0)/(2.0_R8*geo%cvBb(iCv,3)**2)        
            sigin(iCv,1,is,k) = sigin(iCv,0,is,k)
          enddo
         endif
        enddo
      enddo

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tqin

      end
