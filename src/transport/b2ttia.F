*-----------------------------------------------------------------------
*.specification

      subroutine b2ttia (nx, ny, ns, ti, ne2, lnlam, tauia)           !srv 16.06.08 { 01.07.09
      use b2mod_types
      use b2mod_constants
      use b2mod_rates
      use b2mod_b2cmpa
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     &   ti(-1:nx,-1:ny), ne2(-1:nx,-1:ny), lnlam(-1:nx,-1:ny)        !srv 01.07.09
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   tauia(-1:nx,-1:ny,0:ns-1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TTIA computes the classical Braginskii time of the ion-atom collision.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the classical Braginskii time of the ion-atom collision.
*     This coefficient is differed from Balescu coefficient by a factor sqrt(2).
*     - tauia: the classical Braginskii time of the ion-atom collision.
*
*  4. references
*
*     /1/ R. Balescu: Transport Processes in Plasmas, Vol. 1.
*     North Holland, Amsterdam, 1988.
**
*-----------------------------------------------------------------------
*.declarations

*   ..global named constants
      
*   ..local variables
      integer ncall, n2, ix, iy, is
      real (kind=R8) ::
     *   fki, eta, ctaup, t0, t1
      parameter (fki=2.253e0_R8, eta=1.357e0_R8)
      save ncall
*   ..procedures
      intrinsic sqrt
      external subini, subend, xertst, sfill
      external b2xvsg
*   ..initialization
      data ncall/0/

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
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2ttia')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test state
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
       call b2xvsg (n2, ne2, 1, 'ne2', '.gt.')
      endif
*   ..compute ctaup
*     Following Braginskii                                           !srv 11.07.99
*     The mathematical expression is:
*       t0 = (3/(4*sqrt(pi)))*(4*pi*eps0/qe**2)**2*ev*sqrt(ev)       !srv delete sqrt(2)
*       ctaup = t0*sqrt(mp)/lnlam
*     The numerical value in SI units is: ctaup = 1.4743e13/lnlam.
*     The following calculation of t0 is designed to avoid overflow
*     even for IEEE single precision arithmetic.
      t1 = (4*pi*eps0/qe)/qe
      t0 = (3/(4*sqrt(pi)))*t1*(t1*ev*sqrt(ev))                      !srv delete sqrt(2)
      ctaup = t0*sqrt(mp) !/lnlam           !srv 01.07.09    it is taken below

* ..compute tauia
      do is = 0, ns-1
       if (.not.is_neutral(is)) then                            
        do iy = -1, ny
         do ix = -1, nx
          t0 = ti(ix,iy)/ev
          tauia(ix,iy,is) = ctaup/lnlam(ix,iy)*sqrt(am(is))*t0*sqrt(t0)/  !srv 01.07.09
     &      (rz2(ix,iy,is)*ne2(ix,iy))
         enddo
        enddo
       else
        tauia(:,:,is) = 0.0e0_R8
       endif
      enddo

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2ttia

      end                                                            !srv 16.06.08 }
