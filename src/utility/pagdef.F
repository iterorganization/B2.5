      subroutine pagdef ()
      use b2mod_types
      implicit none
#ifndef PLPLOT
*     -----------------------------------------------------------------
*     PAGDEF sets GKS and NCAR properties to default values. It should
*     normally be called at the start of each new plot.
*     -----------------------------------------------------------------
      external gspli, gspmi, gstxi, gsfai, gsln, gslwsc_d, gsplci,
     &  gsmk, gsmksc_d, gspmci, gstxfp, gschxp_d, gschsp_d, gstxci,
     &  gsfais, gsfasi, gsfaci, gschh_d
      external reset, set, cprset, cpseti, cpsetr_d, cpsetc
*     -----------------------------------------------------------------
*   ..set GKS parameters to default
*     (four bundle indices)
      call gspli (1)
      call gspmi (1)
      call gstxi (1)
      call gsfai (1)
*     (thirteen individual attributes)
      call gsln (1)
      call gslwsc_d (1.0_R8)
      call gsplci (1)
      call gsmk (1)
      call gsmksc_d (1.0_R8)
      call gspmci (1)
      call gstxfp (1, 0)
      call gschxp_d (1.0_R8)
      call gschsp_d (0.0_R8)
      call gstxci (1)
      call gsfais (1)
      call gsfasi (1)
      call gsfaci (1)
*     (revised default character height)
      call gschh_d (1.4e-2_R8)
*   ..set NCAR properties to default
*     (I am not sure whether subroutine "reset" does anything)
      call reset ()
*   ..specify default window and viewport
*     (The following call to "set" implies GKS calls to gsvp, gswn,
*     and gselnt for GKS transformation number 1)
      call set_d (0.0_R8, 1.0_R8, 0.0_R8, 1.0_R8, 0.0_R8, 
     & 1.0_R8, 0.0_R8, 1.0_R8, 1)
*   ..initialise parameters for CONPACK
*     (first set everything to default)
      call cprset ()
*     (cls=0: calling program must specify contour levels)
      call cpseti ('CLS', 0)
*     (hls: appearance of high/low labels)
      call cpsetr_d ('HLS', 0.025_R8)
*     (ilt=' ': no information label beneath the plot)
      call cpsetc ('ILT', ' ')
*     (llo, llp, lls: appearance of line labels)
      call cpseti ('LLO', 1)
      call cpseti ('LLP', 3)
      call cpsetr_d ('LLS', 0.016_R8)
*     (nel, net, neu, nls, nlz, nof, nsd: appearance of real*8 numbers)
      call cpseti ('NEL', 2)
      call cpseti ('NET', 0)
      call cpseti ('NEU', 0)
      call cpseti ('NLS', 0)
      call cpseti ('NLZ', 1)
      call cpseti ('NOF', 0)
      call cpseti ('NSD', 3)
*     (set=0: conpack will not call routine set)
      call cpseti ('SET', 0)
#endif
      return
*     -----------------------------------------------------------------
      end subroutine pagdef
