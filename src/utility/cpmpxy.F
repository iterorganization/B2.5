#ifndef GRREAL
#define GRREAL 4
#endif
      subroutine cpmpxy (imap, xinp, yinp, xotp, yotp)
      use b2mod_types
      implicit none
      integer imap
      real (kind=R4) xinp, yinp, xotp, yotp
*     ------------------------------------------------------------------
*     CPMPXY is a substitute for a routine of the same name in the
*     CONPACK sub-package of NCARG graphics. It defines a mapping from
*     rectangular array index coordinates in (1:*,1:*) to the plotting
*     window.
*     ------------------------------------------------------------------
      external gnmap2, gnmap3
*     ------------------------------------------------------------------
*     This routine interfaces to gnmap2 or gnmap3 depending on the
*     value of (imap). Also, the location (xinp,yinp) is shifted to
*     (xinp-1,yinp-1) to account for a different array indexing style.
*     ------------------------------------------------------------------
      if (imap.eq.2) then
       call gnmap2 (xinp-1, yinp-1, xotp, yotp)
      else if (imap.eq.3) then
       call gnmap3 (xinp-1, yinp-1, xotp, yotp)
      else
       xotp = xinp
       yotp = yinp
      endif
      return
*     ------------------------------------------------------------------
      end
