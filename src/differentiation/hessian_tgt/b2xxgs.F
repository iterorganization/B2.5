*-----------------------------------------------------------------------
*.specification

      subroutine b2xxgs (nx, ny, val, fun, iloc)
      use b2mod_types
      use b2mod_b2cmgs
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, iloc
      real (kind=R8) ::
     *   val
*   ..input/output arguments
      real (kind=R8) ::
     *   fun(-1:nx,-1:ny)
*   ..common blocks

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2XXGS sets a mesh function to a constant value on the singular
*     regions of the B2 geometry.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ireg, ix, iy
*   ..procedures
      intrinsic mod
      external xertst

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
*$$$  call subini ('b2xxgs')
*   ..test nx, ny, iloc
      call xertst (0.le.nx.and.0.le.ny,
     &  'b2xxgs--faulty argument nx, ny')
      call xertst (0.le.iloc.and.iloc.lt.4,
     &  'b2xxgs--faulty argument iloc')

* ..reset specified regions
*   ..loop over regions
      do ireg = 0, nlreg-1
       if (iloc.eq.nlloc(ireg)) then
*     ..test region boundaries
        call xertst (
     &    -1+mod(iloc,2).le.nlxlo(ireg).and.nlxhi(ireg).le.nx.and.
     &    -1+iloc/2.le.nlylo(ireg).and.nlyhi(ireg).le.ny,
     &    'b2xxgs--faulty region bounds in /b2cmgs/')
*     ..reset fun on region
        do iy = nlylo(ireg), nlyhi(ireg)
         do ix = nlxlo(ireg), nlxhi(ireg)
          fun(ix,iy) = val
         enddo
        enddo
       endif
      enddo

* ..return
*$$$  call subend ()
      return

*-----------------------------------------------------------------------
*.end b2xxgs

      end subroutine b2xxgs
