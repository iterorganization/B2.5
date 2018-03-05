*-----------------------------------------------------------------------
*.specification

      subroutine b2treq (nx, ny, vol, te, ti, ne, ni, ne2m, lnlam, ceqp)   !xpb !srv 01.07.09
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_external
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny
      real (kind=R8) ::
     *  vol(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1), 
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny),
     &  ne(-1:nx,-1:ny), ne2m(-1:nx,-1:ny), lnlam(-1:nx,-1:ny)        !xpb !srv 01.07.09
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  ceqp(-1:nx,-1:ny)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TREQ computes the cell-integrated energy equipartition
*     multiplier.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the electron-atom energy equipartition
*     multiplier ceqp, defined so that on each cell the signed
*     equipartition energy transfer from atoms to electrons is
*     (ceqp*ne*ni/te**1.5)*(te-ti) .
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy
      integer, save :: icase = 0
      real (kind=R8) ::
     *  cfeqp, phm0
      parameter (cfeqp=4.8e-15_R8)
      save ncall, phm0
*     lnlam is the log of the plasma parameter.
*     cfeqp is a numerical coefficient in the energy transfer rate.
*   ..procedures
      intrinsic sqrt
      external subini, subend, xertst, ipgetr
      external b2xvsg
*   ..initialisation
      data ncall/0/, phm0/1.0_R8/

*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2treq')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..read input data
      if (ncall.eq.0) then
       call ipgetr ('b2treq_phm0', phm0)
       call ipgeti ('b2tlnl_ei', icase)
       call xertst (0.eq.icase .or. 3.eq.icase, 
     &              'faulty input b2tlnl_ei') 
       call xertst (0.le.phm0, 'faulty input phm0') 
      endif
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test vol, state
       call b2xvsg (n2, vol, 1, 'vol', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2*2, ni, 1, 'ni', '.gt.')
       call b2xvsg (n2, ne2m, 1, 'ne2m', '.gt.')
      endif

      call b2tlnl (nx, ny, icase, te, ti, ne, lnlam)                     !srv 20.09.11
*   ..compute ceqp
      do iy = -1, ny
       do ix = -1, nx
        ceqp(ix,iy) = phm0*(cfeqp*lnlam(ix,iy)*ev*sqrt(ev))*vol(ix,iy)*  !srv 01.07.09
     &    (ne2m(ix,iy)/(ni(ix,iy,0)+ni_ext(ix,iy,0)))
       enddo
      enddo

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2treq

      end subroutine b2treq

