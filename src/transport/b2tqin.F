*-----------------------------------------------------------------------
*.specification

      subroutine b2tqin (nx, ny, ns, nscx, iscx, bb, ti, na, sigin)
      use b2mod_types
      use b2mod_rates
      use b2mod_constants
      use b2mod_b2cmpa
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, nscx
      integer iscx(0:nscx-1)
      real (kind=R8) ::
     *  bb(-1:nx,-1:ny,0:3), ti(-1:nx,-1:ny),                    !xpb
     &  na(-1:nx,-1:ny,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *	sigin(-1:nx,-1:ny,0:1,0:ns-1,0:nscx-1)                   !sv 01.07.99 !xpb
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
      integer ncall, n2, ix, iy, is, is0, k
      real (kind=R8) ::
     *   phm0, cutlo, cutll, rcx, t0
      save ncall, phm0
*   ..procedures
      intrinsic sqrt, log, exp, max, dim
      real (kind=R8) ::
     *   machsfr, expu
      external subini, subend, xertst, sfill, ipgetr, machsfr
      external b2xvsg
      expu(t0) = dim(exp(max(t0,cutll)),cutlo)
*     (expu is the exponential function with underflow set to 0)
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
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..set cutll, cutlo
      cutll = log(machsfr())
      cutlo = max(machsfr(),exp(cutll))
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test bb
       call b2xvsg (n2, bb(-1,-1,3), 1, 'bb3', '.gt.')
*    ..test state
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
       call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
      endif
    
* ..compute classical sigin
*   ..compute sigin.
      call sfill (n2*2*ns*nscx, 0.0_R8, sigin, 1)
*   ..we use the charge exchange rate as computed from b2frates
*   ..this is only nonzero for hydrogen neutral on hydrogen ion CX
      do k = 0, nscx-1
        is0 = iscx(k)
        do is = 0, ns-1
         if (.not.is_neutral(is).and.nint(zn(is)).eq.1) then
          do iy = -1, ny
           do ix = -1, nx
            t0 = log(ti(ix,iy)/(am(is0)*ev))
            rcx = phm0*
     &       expu(rlcx(ix,iy,0,is,k)+rlcx(ix,iy,1,is,k)*t0)
! From 5.2: Using the diffusivity of room temperature hydrogen gas
            sigin(ix,iy,0,is,k) = (na(ix,iy,is)*am(is)*mp)      !srv 01.07.99
     &       *(na(ix,iy,is0)*3.2e-15_R8)
     &       *sqrt(ti(ix,iy)/ev/0.026_R8)/(2.0_R8*bb(ix,iy,3)**2)
! From 5.0/5.1: Using the CX rate
!            sigin(ix,iy,0,is,k) = (na(ix,iy,is)*am(is)*mp)
!     &       *rcx*na(ix,iy,is0)/(2.0_R8*bb(ix,iy,3)**2)        
            sigin(ix,iy,1,is,k) = sigin(ix,iy,0,is,k)
           enddo
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
