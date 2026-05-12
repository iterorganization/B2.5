*-----------------------------------------------------------------------
*.specification

      subroutine b2rurc_mol (nget)
      use b2mod_types
      use b2mod_math
      use b2mod_b2cmrc_mol
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nget
cc      character*10 version !! re-introduce versions once initial development phase is over. 
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2RURC_MOL reads a table of molecular rate data from an un*formatted
*     file.
*
*  3. description (see also routine b2cdca)
*
*     This routine reads from unit (nget) a table of molecular rate data.
*     Unit nget must be connected to an un*formatted file.
*     See the code to discern the required format of the input file.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nget - integer, input.
*     nget specifies a unit number for un*formatted input.
*
*     (parameters in common:)
*
*     /b2cmrc_mol/ - output.
*     /b2cmrc_mol/ contains a table of molecular rate coefficients.
*
*      rtnt - integer constant.
*     rtnt specifies the number of intervals in the discretisation
*     of the temperature range. The table will contain data for
*     (rtnt+1) values of the temperature variable.
*
*      rtnn - integer constant.
*     rtnn specifies the number of intervals in the discretisation
*     of the density range. The table will contain data for (rtnn+1)
*     values of the density variable.
*
*      rtnms - integer.
*     rtnms specifies the total number of molecular species. 
*      rtzmin, rtzmax, rtzn - (0:rtnms-1) real*8 array.
*
*      rtt - (0:rtnt) real*8 array.
*     rtt(0:rtnt) specifies the sequence of values of the temperature
*     variable, expressed in ev, for which the rate coefficients are
*     evaluated. rtt(0:rtnt) will be positive and increasing.
*
*      rtn - (0:rtnn) real*8 array.
*     rtn(0:rtnn) specifies the sequence of values of the density
*     variable, expressed in m**-3, for which the rate coefficients
*     are evaluated. rtn(0:rtnn) will be positive and increasing.
*
*      rtlt - (0:rtnt) real*8 array.
*     rtlt(0:rtnt) holds log(rtt(0:rtnt)).
*
*      rtln - (0:rtnn) real*8 array.
*     rtln(0:rtnn) holds log(rtn(0:rtnn)).
*
*      rtl_nd - (0:rtnt,0:rtnn,0:rtnms-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtnms-1), rtl_nd(i,j,k) specifies
*     logu(rt_nd), where rt_nd is the non-dissociative ionisation rate coefficient at
*     electron temperature rtt(i), electron density rtn(j).
*     rt_d has dimension m**3/s. The function logu is
*     the same as (log) except that logu(0.0) is equal to some huge
*     negative number.
*
*
*  8. further comments
*
*     See routine b2cdca for explanation of 'un*formatted file'.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer idum(0:9)
      real (kind=R8) ::
     *   lnzero
      parameter (lnzero=-1.0e30_R8)
*   ..procedures
      real (kind=R8) ::
     *   smax
      external xertst, smax,
     *   cfrure, cfrure_dep, cfruin, cfruin_opt, cfruin_dep

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2rurc_mol')
*   ..test nget
      call xertst (1.le.nget, 'faulty argument nget')
*   ..set cutll
      if (.not.b2mod_math_initialised) call b2mod_math_init

* ..read rate coefficients

*   ..read dimensions
      call cfruin (nget, 3, idum(0:2), 
     .    'rtnt,rtnn,rtnms') 
      call cfruin (nget, 2, idum(3:4), 
     .    'rtntr,rtn_mcoll')
      call alloc_b2mod_b2cmrc_mol(idum(0), idum(1), idum(2),
     .     idum(3), idum(4))
*   ..read species data

*   ..read abscissae
cc values
      call cfrure (nget, rtnt+1, rtt, 'rtt')
      call cfrure (nget, rtnn+1, rtn, 'rtn')
      call cfrure (nget, rtntr+1, rttr, 'rttr')
      call cfrure (nget, rtn_mcoll+1, rt_mcoll, 'rt_mcoll')
cc logs of values
      call cfrure (nget, rtnt+1, rtlt, 'rtlt')
      call cfrure (nget, rtnn+1, rtln, 'rtln')
      call cfrure (nget, rtntr+1, rtltr, 'rtltr')
      call cfrure (nget, rtn_mcoll+1, rtl_mcoll, 'rtl_mcoll')
*   ..read rate coefficients
      call cfrure (nget, (rtnt+1)*(rtnn+1)*rtnms, rtl_nd, 'rtl_nd')
      call cfrure (nget, (rtnt+1)*(rtnn+1)*rtnms, rtl_d,  'rtl_d')
      call cfrure (nget, (rtnt+1)*(rtnn+1)*rtnms, rtl_di1,'rtl_di1')
      call cfrure (nget, (rtnt+1)*(rtnn+1)*rtnms, rtl_ic, 'rtl_ic')
      call cfrure (nget, (rtnt+1)*(rtnn+1)*rtnms, rtl_dr, 'rtl_dr')
      call cfrure (nget, (rtnt+1)*(rtnn+1)*rtnms, rtl_de, 'rtl_de')
      call cfrure (nget, (rtnt+1)*(rtnn+1)*rtnms, rtl_di2,'rtl_di2')
      call cfrure (nget, (rtnt+1)*(rtntr+1)*(rtn_mcoll+1)*
     .            rtnms, rtl_ec_m, 'rtl_ec_m')
      call cfrure (nget, (rtnt+1)*(rtntr+1)*
     .            rtnms, rtl_ec_E, 'rtl_ec_E')
* ..test rate coefficients
 !! TODO

      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2rurc
      end subroutine b2rurc_mol
