      subroutine pagend
      use b2mod_types
      implicit none
*     -----------------------------------------------------------------
*     PAGEND terminates a plot.
*     This version is for NCARGKS.
*     -----------------------------------------------------------------
#ifndef PLPLOT
      external frame
#endif
*     -----------------------------------------------------------------
*     (The following call to "frame" implies a GKS calls to "gclrwk"
*     for each open workstation.)
#ifndef PLPLOT
      call frame ()
#endif
      return
*     -----------------------------------------------------------------
      end subroutine pagend
