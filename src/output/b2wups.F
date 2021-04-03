*-----------------------------------------------------------------------
*.specification

      subroutine b2wups (nput, nx, ny, ns,
     &  ne, na, ua, uadia, te, ti, tn, po, kt, zt,
     &  fna, fhe, fhi, fhn, fch, fch_32, fch_52, fch_p, kinrgy,
     &  fna_mdf, fhe_mdf, fhi_mdf, fna_fcor, fna_nodrift, fna_he,
     &  fnaPSch, fhePSch, fhiPSch, fna_eir, fne_eir, fhe_eir, fhi_eir,
     &  fna_32, fna_52, fni_32, fni_52, fne_32, fne_52,
     &  fchdia, fchin, fchvispar, fchvisper, fchvisq, fchinert,
     &  vaecrb, vadia, wadia, veecrb, vedia, floe_noc, floi_noc)
      use b2mod_types
      use b2mod_time
      use b2mod_version
      implicit none
*   ..input arguments (unchanged on exit)
      integer nput, nx, ny, ns
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  ne(-1:nx,-1:ny), uadia(-1:nx,-1:ny,0:1,0:1,0:ns-1),                  !xpb
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), tn(-1:nx,-1:ny),
     &  po(-1:nx,-1:ny), kt(-1:nx,-1:ny), zt(-1:nx,-1:ny)
      real (kind=R8) ::
     *  fna(-1:nx,-1:ny,0:1,0:1,0:ns-1), fhe(-1:nx,-1:ny,0:1,0:1),
     &  fhi(-1:nx,-1:ny,0:1,0:1), fhn(-1:nx,-1:ny,0:1,0:1), 
     &  fch(-1:nx,-1:ny,0:1,0:1),
     1  fch_32(-1:nx,-1:ny,0:1,0:1), fch_52(-1:nx,-1:ny,0:1,0:1),
     2  fch_p(-1:nx,-1:ny,0:1,0:1), kinrgy(-1:nx,-1:ny,0:ns-1),          !srv 10.06.08
     &  fna_mdf(-1:nx,-1:ny,0:1,0:1,0:ns-1), 
     &  fhe_mdf(-1:nx,-1:ny,0:1,0:1),
     &  fhi_mdf(-1:nx,-1:ny,0:1,0:1), 
     &  fna_fcor(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     &  fna_nodrift(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     &  fna_he(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     &  fnaPSch(-1:nx,-1:ny,0:1,0:1,0:ns-1), 
     &  fhePSch(-1:nx,-1:ny,0:1,0:1),
     &  fhiPSch(-1:nx,-1:ny,0:1,0:1), 
     &  fna_eir(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     &  fne_eir(-1:nx,-1:ny,0:1,0:1), fhe_eir(-1:nx,-1:ny,0:1,0:1),
     &  fhi_eir(-1:nx,-1:ny,0:1,0:1),
     &  fna_32(-1:nx,-1:ny,0:1,0:1,0:ns-1), 
     &  fna_52(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     &  fni_32(-1:nx,-1:ny,0:1,0:1), fni_52(-1:nx,-1:ny,0:1,0:1),
     &  fne_32(-1:nx,-1:ny,0:1,0:1), fne_52(-1:nx,-1:ny,0:1,0:1),
     &  fchdia(-1:nx,-1:ny,0:1,0:1), fchin(-1:nx,-1:ny,0:1,0:1),
     &  fchvispar(-1:nx,-1:ny,0:1,0:1), fchvisper(-1:nx,-1:ny,0:1,0:1),
     &  fchvisq(-1:nx,-1:ny,0:1,0:1), fchinert(-1:nx,-1:ny,0:1,0:1),
     &  vaecrb(-1:nx,-1:ny,0:1,0:1,0:ns-1), 
     &  vadia(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     &  wadia(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     &  veecrb(-1:nx,-1:ny,0:1,0:1), vedia(-1:nx,-1:ny,0:1,0:1),
     &  floe_noc(-1:nx,-1:ny,0:1,0:1), floi_noc(-1:nx,-1:ny,0:1,0:1)
*   ..output arguments (unspecified on entry)
*     (none)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WUPS writes the state of the plasma to an un*formatted file.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine writes to unit (nput) the state of the plasma:
*     fields na, ua, te, ti, po, and fluxes fna, fhe, fhi, fch.
*     Unit nput must be connected to an un*formatted file. See the
*     code to discern the format of the output file.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nput - integer, input.
*     nput specifies a unit number for un*formatted output.
*
*     nx, ny, ns - integers, input.
*     See the calling routine for description.
*
*     ne, na, ua, uadia, te, ti, po, fna, fhe, fhi, fch - real*8 arrays, input.
*     See the calling routine for description.
*
*
*  8. further comments
*
*     See routine b2cdca for explanation of 'un*formatted file'.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer n2
      real (kind=R8) :: time(1)
*   ..procedures
      external subini, subend, xertst, cfwure
      external b2xvps, b2xvsg

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wups')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nput
!WG      call xertst (1.le.nput, 'faulty argument nput')
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')

* ..save state and fluxes
*   ..save plasma state
      call cfwure (nput, n2*ns, na, 'na')
      call cfwure (nput, n2,    ne, 'ne')
      call cfwure (nput, n2*ns, ua, 'ua')
      call cfwure (nput, n2*4*ns, uadia, 'uadia')
      call cfwure (nput, n2, te, 'te')
      call cfwure (nput, n2, ti, 'ti')
      call cfwure (nput, n2, tn, 'tn')
      call cfwure (nput, n2, po, 'po')
      call cfwure (nput, n2, kt, 'kt')
      call cfwure (nput, n2, zt, 'zt')
*   ..save fluxes
      call cfwure (nput, n2*4*ns, fna, 'fna')
      call cfwure (nput, n2*4, fhe, 'fhe')
      call cfwure (nput, n2*4, fhi, 'fhi')
      call cfwure (nput, n2*4, fhn, 'fhn')
      call cfwure (nput, n2*4, fch, 'fch')
      call cfwure (nput, n2*4, fch_32, 'fch_32')
      call cfwure (nput, n2*4, fch_52, 'fch_52')
      call cfwure (nput, n2*ns, kinrgy, 'kinrgy')
*   ..save time
      time(1)=tim
      call cfwure (nput, 1, time, 'time')
*   ..save parallel current
      call cfwure (nput, n2*4, fch_p, 'fch_p')
*   ..saving various fluxes and velocities needed by advanced SOLPS-ITER features
      call cfwure (nput, n2*4*ns, fna_mdf, 'fna_mdf')
      call cfwure (nput, n2*4, fhe_mdf, 'fhe_mdf')
      call cfwure (nput, n2*4, fhi_mdf, 'fhi_mdf')
      call cfwure (nput, n2*4*ns, fna_fcor, 'fna_fcor')
      call cfwure (nput, n2*4*ns, fna_nodrift, 'fna_nodrift')
      call cfwure (nput, n2*4*ns, fna_he, 'fna_he')
      call cfwure (nput, n2*4*ns, fnaPSch, 'fnaPSch')
      call cfwure (nput, n2*4, fhePSch, 'fhePSch')
      call cfwure (nput, n2*4, fhiPSch, 'fhiPSch')
      call cfwure (nput, n2*4*ns, fna_eir, 'fna_eir')
      call cfwure (nput, n2*4, fne_eir, 'fne_eir')
      call cfwure (nput, n2*4, fhe_eir, 'fhe_eir')
      call cfwure (nput, n2*4, fhi_eir, 'fhi_eir')
      call cfwure (nput, n2*4*ns, fna_32, 'fna_32')
      call cfwure (nput, n2*4*ns, fna_52, 'fna_52')
      call cfwure (nput, n2*4, fni_32, 'fni_32')
      call cfwure (nput, n2*4, fni_52, 'fni_52')
      call cfwure (nput, n2*4, fne_32, 'fne_32')
      call cfwure (nput, n2*4, fne_52, 'fne_52')
      call cfwure (nput, n2*4, fchdia, 'fchdia')
      call cfwure (nput, n2*4, fchin, 'fchin')
      call cfwure (nput, n2*4, fchvispar, 'fchvispar')
      call cfwure (nput, n2*4, fchvisper, 'fchvisper')
      call cfwure (nput, n2*4, fchvisq, 'fchvisq')
      call cfwure (nput, n2*4, fchinert, 'fchinert')
      call cfwure (nput, n2*4*ns, vaecrb, 'vaecrb')
      call cfwure (nput, n2*4*ns, vadia, 'vadia')
      call cfwure (nput, n2*4*ns, wadia, 'wadia')
      call cfwure (nput, n2*4, veecrb, 'veecrb')
      call cfwure (nput, n2*4, vedia, 'vedia')
      call cfwure (nput, n2*4, floe_noc, 'floe_noc')
      call cfwure (nput, n2*4, floi_noc, 'floi_noc')

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2wups

      end subroutine b2wups
