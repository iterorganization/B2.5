      subroutine gncolr (acrx, acry, nacra, iaia, igia, naia)
      use b2mod_types
      implicit none
      integer nacra, naia
      real (kind=R4) acrx(0:nacra-1), acry(0:nacra-1)
      integer iaia(0:naia-1), igia(0:naia-1)
*     ------------------------------------------------------------------
*     GNCOLR is for use by the routine ARSCAM in the AREAS utility of
*     NCAR graphics. It defines the filling of a polygon.
*     ------------------------------------------------------------------
      logical dofill
      integer i
      external lbfill
*     ------------------------------------------------------------------
*     The routine follows the "colram" example in the NCAR-3 manual,
*     Contours page 3-199.
*     ------------------------------------------------------------------
*     (Do not fill if any area identifier is negative)
      dofill = .true.
      do i = 0, naia-1
       dofill = dofill.and.0.le.iaia(i)
      enddo
*     If (dofill), search for group identifier 3 (which is associated
*     with the contour lines) and fill according to the corresponding
*     area identifier. Note that group 3 is not necessarily present.
      if (dofill) then
       i = 0
    1  if (i.lt.naia) then
        if (igia(i).ne.3) then
         i = i+1
         goto 1
        endif
       endif
       if (i.lt.naia) then
        call lbfill (acrx, acry, nacra, iaia(i))
       endif
      endif
      return
*     ------------------------------------------------------------------
      end subroutine gncolr
