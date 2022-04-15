*-----------------------------------------------------------------------
*.specification

      subroutine b2nxfv (nx, ny, isb, ismain, hz, cvsbhz_cl, flubv)    !srv 19.01.05 01.07.05 11.09.09 { 16.10.17 19.02.21
      use b2mod_types
      use b2mod_indirect
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, isb, ismain                                      !srv 19.02.21
      real (kind=R8) ::
     *   hz(-1:nx,-1:ny), cvsbhz_cl(-1:nx,-1:ny,0:1,0:1)               !srv 16.10.17
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   flubv(-1:nx,-1:ny,0:1,0:1)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2NXFV calculates the convective term due to new form of ion  
*     viscosity term in the momentum correction equation.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, ix, iy
      real (kind=R8) ::
     *   ft                                                            !srv 16.10.17
      save ncall
      !$OMP THREADPRIVATE(ncall)
      external xertst, ipgeti, ipgetr, b2scopy
      external b2xvff, b2xvfx, b2xvfv, b2xvsg
*   ..initialization
      data ncall/0/

*-----------------------------------------------------------------------
*.documentation-internal
*
*     Description of some local variables:
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2nxfv')
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
      endif

* ..initialize flubv to zero
      flubv = 0.0_R8

* ..compute convection term
      ft = 4.0_R8/3.0_R8
*   ..compute flubvx
!      if(nint(zn(isb)).eq.1 .and. nint(zamax(isb)).eq.1) then         !srv 26.10.08
      if(isb.eq.ismain) then                                          !srv 19.02.21
       do iy = -1, ny
        do ix = -1, nx
         if(leftix(ix,iy).ne.-2) then
          flubv(ix,iy,0,0) = ft*cvsbhz_cl(ix,iy,0,0)*                 !srv 16.10.17
     &      (log(hz(ix,iy))-log(hz(leftix(ix,iy),leftiy(ix,iy))))
         endif
        enddo
       enddo
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2nxfv

      end subroutine b2nxfv                                           !srv 19.01.05 01.07.05 11.09.09 }
