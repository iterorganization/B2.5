*-----------------------------------------------------------------------
*.specification

      subroutine b2nxfv (nx, ny, isb, hz, cvsb, flubv)               !srv 19.01.05 01.07.05 11.09.09 {
      use b2mod_types
      use b2mod_indirect
      use b2mod_b2cmpa
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, isb
      real (kind=R8) ::
     *   hz(-1:nx,-1:ny), cvsb(-1:nx,-1:ny,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   flubv(-1:nx,-1:ny,0:1)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2NXFV calculates the convective term due to new form of ion  
*     viscosity term in the momentum correction  equation.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy
      save ncall
      external subini, subend, xertst, ipgeti, ipgetr, scopy
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
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
      endif

* ..compute convection term
*   ..compute flubvx
      if(nint(zn(isb)).eq.1 .and. nint(zamax(isb)).eq.1) then         !srv 26.10.08
       do iy = -1, ny
        do ix = -1, nx
	   if(leftix(ix,iy).ne.-2) then
          flubv(ix,iy,0) = cvsb(ix,iy,0)*
     &      (log(hz(ix,iy))-log(hz(leftix(ix,iy),leftiy(ix,iy))))
	   else
          flubv(ix,iy,0) = 0.0e0_R8
	   endif
        enddo
       enddo

	 flubv(:,:,1) = 0.0_R8
	else
	 flubv = 0.0_R8
	endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2nxfv

      end                                                            !srv 19.01.05 01.07.05 11.09.09 }
