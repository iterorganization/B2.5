*-----------------------------------------------------------------------
*.specification

      subroutine b2wurc_mol (nput, version)
      use b2mod_types
      use b2mod_math
      use b2mod_b2cmrc_mol
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nput
      character*10 version
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WURC_MOL writes a table of molecular rate coefficients to an
*     un*formatted file.
*
*
*  3. description 
*
*     This routine writes to unit (nput) a table of molecular rate
*     coefficients.
*     Unit nput must be connected to an un*formatted file.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nput - integer, input.
*     nput specifies a unit number for un*formatted output.
*
*     (parameters in common:)
*
*
*  8. further comments
*
*     See routine b2cdca for explanation of 'un*formatted file'.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      logical lt0
      integer k, idum(0:9)
      real (kind=R8) ::
     *   t0
*   ..procedures
      real (kind=R8) ::
     *   smax
      external xertst, cfwure, cfwuin
      external smax

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wurc_mol')
*   ..set cutll
      if (.not.b2mod_math_initialised) call b2mod_math_init
*   ..test nput
      call xertst (1.le.nput, 'faulty argument nput')

* ..test rate coefficients
 !! TODO

* ..save the table of rate coefficients
*   ..save dimensions
      idum(0) = rtnt
      idum(1) = rtnn
      idum(2) = rtnms
      idum(3) = rtntr
      idum(4) = rtn_mcoll

      call cfwuin (nput, 3, idum(0:2), 
     .    'rtnt,rtnn,rtnms')
      call cfwuin (nput, 2, idum(3:4), 
     .    'rtntr,rtn_mcoll')

*   ..save abscissae
      call cfwure (nput, rtnt+1, rtt, 'rtt')
      call cfwure (nput, rtnn+1, rtn, 'rtn')
      call cfwure (nput, rtntr+1, rttr, 'rttr')
      call cfwure (nput, rtn_mcoll+1, rt_mcoll, 'rt_mcoll')
      call cfwure (nput, rtnt+1, rtlt, 'rtlt')
      call cfwure (nput, rtnn+1, rtln, 'rtln')
      call cfwure (nput, rtntr+1, rtltr, 'rtltr')
      call cfwure (nput, rtn_mcoll+1, rtl_mcoll, 'rtl_mcoll')

*   ..save (log of) rate coefficients
      call cfwure (nput, (rtnt+1)*(rtnn+1)*rtnms, rtl_nd, 'rtl_nd')
      call cfwure (nput, (rtnt+1)*(rtnn+1)*rtnms, rtl_d, 'rtl_d')
      call cfwure (nput, (rtnt+1)*(rtnn+1)*rtnms, rtl_di1, 'rtl_di1')
      call cfwure (nput, (rtnt+1)*(rtnn+1)*rtnms, rtl_ic, 'rtl_ic')
      call cfwure (nput, (rtnt+1)*(rtnn+1)*rtnms, rtl_dr, 'rtl_dr')
      call cfwure (nput, (rtnt+1)*(rtnn+1)*rtnms, rtl_de, 'rtl_de')
      call cfwure (nput, (rtnt+1)*(rtnn+1)*rtnms, rtl_di2, 'rtl_di2')
      call cfwure (nput, (rtnt+1)*(rtntr+1)*(rtn_mcoll+1)*
     .     rtnms, rtl_ec_m, 'rtl_ec_m')
      call cfwure (nput, (rtnt+1)*(rtntr+1)*
     .     rtnms, rtl_ec_E, 'rtl_ec_E')

* ..return
      flush (nput)
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2wurc

      end subroutine b2wurc_mol
