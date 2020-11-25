*-----------------------------------------------------------------------
*.specification

      subroutine b2nxfm (nx, ny, mb, fnb, flub)
      use b2mod_types
      use b2mod_indirect
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny
      real (kind=R8) ::
     *   mb, fnb(-1:nx,-1:ny,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   flub(-1:nx,-1:ny,0:1)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2NXFM computes the flow terms in the parallel momentum equation.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine evaluates the coefficients flubx and fluby that
*     occur in the expression for the transport of parallel momentum
*     for one selected species; symbolically:
*       fmbx = flubx*ub-cvsbx*(d/dx)ub
*       fmby = fluby*ub-cvsby*(d/dy)ub
*     cvsbx and cvsby are viscosities (not evaluated in this routine),
*     and the expressions for flubx and fluby are:
*       flubx = mb*fnbx
*       fluby = mb*fnby
*     where mb is the species mass and fnbx and fnby are particle
*     fluxes. Some interpolation between cell centers, faces, and
*     corners is involved, and one is referred to the code for the
*     precise expressions.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, ix, iy
      save ncall
*   ..procedures
      external subini, subend, xertst, sfill
      external b2xvff
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2nxfm')
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..test mb
      call xertst (0.lt.mb, 'faulty argument mb')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test edge values of fnb
       call b2xvff (nx, ny, fnb, 'fnb')
      endif

* ..compute convection factors
*   ..compute flubx
      do iy = -1, ny
       do ix = -1, nx
        if(leftix(ix,iy).eq.-2) then
         flub(ix,iy,0) = 0
        elseif(rightix(ix,iy).eq.nx+1) then
         flub(ix,iy,0) = 0
        else
         flub(ix,iy,0) = (mb/2)*
     1    (fnb(ix,iy,0)+fnb(rightix(ix,iy),rightiy(ix,iy),0))
        endif
       enddo
      enddo
*   ..compute fluby
      do iy = -1, ny
       do ix = -1, nx
        if(bottomiy(ix,iy).eq.-2) then
         flub(ix,iy,1) = 0
        elseif(topiy(ix,iy).eq.ny+1) then
         flub(ix,iy,1) = 0
        else
         flub(ix,iy,1) = (mb/2)*
     1    (fnb(topix(ix,iy),topiy(ix,iy),1)+fnb(ix,iy,1))
        endif
       enddo
      enddo

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2nxfm

      end subroutine b2nxfm
