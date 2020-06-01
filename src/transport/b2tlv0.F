      subroutine b2tlv0 (nCv, nFc, nVx, ns, geo, mpg, na, ua, te, 
     &  ti, rza, cvsa0, cvsa, cvsahz, flv)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_map
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *  na(nCv,0:ns-1), ua(nCv,0:ns-1), te(nCv), ti(nCv), 
     *  rza(nCv,0:ns-1), cvsa0(nFc,0:ns-1)
*   ..input and output arguments
      real (kind=R8) ::
     *  cvsa(nFc,0:1,0:ns-1), cvsahz(nFc,0:1,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  flv(nFc,0:1,0:ns-1)
*   ..common blocks
*     ------------------------------------------------------------------
*     B2TLV0 computes the neutral viscosity flux limit.
*     ------------------------------------------------------------------
*   ..local variables
      integer iFc, is, ncall
      integer, save :: style = 0
      real (kind=R8) ::
     *  alpha, gamma, qcl, qfl, tfl, maxfluxlimit(0:1),
     *  dubf(nFc,0:1), ubv(nVx)
      real (kind=R8), save :: cutlo = 1.0e-60_R8
      save ncall, alpha, gamma
*   ..initialization
      data ncall/0/, alpha/0.0_R8/, gamma/2.0_R8/
*   ..procedures
      intrinsic max, abs, sqrt
      external subini, subend, xertst, ipgetr
*     ------------------------------------------------------------------
*   ..subprogram start-up calls
      call subini ('b2tlv0')
*   ..test nCv, nFc, ns
      call xertst (0.le.nCv .and. 0.le.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..read input on first call
      if (ncall.eq.0) then
        call ipgetr ('b2tlv0_alpha', alpha)
        call xertst (0.le.alpha, 'invalid input parameter alpha')
        call ipgetr ('b2tlv0_gamma', gamma)
        call xertst (0.lt.gamma, 'invalid input parameter gamma')
        call ipgeti ('b2tlv0_style', style)
      endif

      call sfill(nFc*2*ns, 1.0_R8, flv, 1)
      call sfill(2, 1.0_R8, maxfluxlimit, 1)

*   ..nothing done if alpha.eq.0
      if (alpha.ne.0.0_R8) then
       if (style.eq.0) then
        do is = 0, ns-1
         if (is_neutral(is)) then
*    ..compute differences of parallel velocity
          call diff (nCv, nFc, nVx, 0, geo, mpg, ua(1,is), ubv, dubf)
          do iFc = 1, nFc
*    ..poloidal
           qcl = abs(cvsa(iFc,0,is)*dubf(iFc,0))
           qfl = alpha*geo%fcS(iFc)*abs(geo%fcQalf(iFc,0))*
     &          (na(mpg%fcCv(iFc,1),is)*ti(mpg%fcCv(iFc,1))+
     &           na(mpg%fcCv(iFc,2),is)*ti(mpg%fcCv(iFc,2)))/2.0_R8
           tfl = 1.0_R8/(1.0_R8 + abs(qcl/(qfl+cutlo))**gamma)
     &                                                 **(1.0_R8/gamma)
            cvsa(iFc,0,is)   = cvsa(iFc,0,is)*tfl
            cvsahz(iFc,0,is) = cvsahz(iFc,0,is)*tfl
            flv(iFc,0,is) = tfl
            maxfluxlimit(0)=max(maxfluxlimit(0),tfl)
*    ..perpendicular
           qcl = abs(cvsa(iFc,1,is)*dubf(iFc,1))
           qfl = alpha*geo%fcS(iFc)*abs(geo%fcQalf(iFc,1))*
     &          (na(mpg%fcCv(iFc,1),is)*ti(mpg%fcCv(iFc,1))+
     &           na(mpg%fcCv(iFc,2),is)*ti(mpg%fcCv(iFc,2)))/2.0_R8
           tfl = 1.0_R8/(1.0_R8 + abs(qcl/(qfl+cutlo))**gamma)
     &                                                 **(1.0_R8/gamma)
            cvsa(iFc,1,is)   = cvsa(iFc,1,is)*tfl
            cvsahz(iFc,1,is) = cvsahz(iFc,1,is)*tfl
            flv(iFc,1,is) = tfl
            maxfluxlimit(1)=max(maxfluxlimit(1),tfl)
	  enddo
         endif
        enddo
       else
cwdk    flux limit total gradient-driven flow, not individual poloidal
cwdk    or radial contributions, since neutral transport is isotropic
        do is = 0, ns-1
         if (is_neutral(is)) then
*    ..compute differences of parallel velocity
          call diff (nCv, nFc, nVx, 0, geo, mpg, ua(1,is), ubv, dubf)
          do iFc = 1, nFc
*    ..limit total gradient
           qcl = cvsa0(iFc,is)*sqrt(dubf(iFc,0)**2 + dubf(iFc,1)**2)
           qfl = alpha*geo%fcS(iFc)*
     &          (na(mpg%fcCv(iFc,1),is)*ti(mpg%fcCv(iFc,1))+
     &           na(mpg%fcCv(iFc,2),is)*ti(mpg%fcCv(iFc,2)))/2.0_R8
           tfl = 1.0_R8/(1.0_R8 + abs(qcl/(qfl+cutlo))**gamma)
     &                                                 **(1.0_R8/gamma)
*    ..poloidal
            cvsa(iFc,0,is)   = cvsa(iFc,0,is)*tfl
            cvsahz(iFc,0,is) = cvsahz(iFc,0,is)*tfl
            flv(iFc,0,is) = tfl
            maxfluxlimit(0)=max(maxfluxlimit(0),tfl)
*    ..radial
            cvsa(iFc,1,is)   = cvsa(iFc,1,is)*tfl
            cvsahz(iFc,1,is) = cvsahz(iFc,1,is)*tfl
            flv(iFc,1,is) = tfl
            maxfluxlimit(1)=max(maxfluxlimit(1),tfl)
	  enddo
         endif
        enddo
       endif
      endif

!      write(*,'(1x,a,1p,2e10.2)') 'b2tlv0: maxfluxlimit = ',maxfluxlimit

*   ..return
      ncall = ncall + 1
      call subend ()
      return
*     ------------------------------------------------------------------
      end
