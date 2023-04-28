      subroutine write_b2fstate_diff (nput, nCv, nFc, ns, st)
      use b2mod_time
      use b2us_plasma_diff
      use b2mod_subsys
      implicit none
      integer, intent(in) :: nput, nCv, nFc, ns
      type(B2STATE_DIFF), intent(inout) :: st
      integer ::  idum(0:2) 
      real (kind=R8) :: time(1)

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('write_b2fstate_diff')
*   ..test ns
      call xertst (1.le.ns, 'faulty argument ns')
      
* ..save state and fluxes

*   ..save plasma state
      call cfwure (nput, nCv*ns,   st % pl % na, 'na')
      call cfwure (nput, nCv,      st % dv % ne, 'ne')
      call cfwure (nput, nCv*ns,   st % pl % ua, 'ua')
      call cfwure (nput, nFc*2*ns, st % dv % uadia, 'uadia') 
      call cfwure (nput, nCv,      st % pl % te, 'te')
      call cfwure (nput, nCv,      st % pl % ti, 'ti')
      call cfwure (nput, nCv,      st % pl % tn, 'tn')
      call cfwure (nput, nCv,      st % pl % po, 'po')
      call cfwure (nput, nCv,      st % pl % kt, 'kt')

*   ..save fluxes
      call cfwure (nput, nFc*2*ns, st % dv % fna, 'fna')
      call cfwure (nput, nFc*2,    st % dv % fhe, 'fhe')
      call cfwure (nput, nFc*2,    st % dv % fhi, 'fhi')
      call cfwure (nput, nFc*2,    st % dv % fhn, 'fhn')
      call cfwure (nput, nFc*2,    st % dv % fch, 'fch')
      call cfwure (nput, nFc*2,    st % dv % fch_32, 'fch_32')
      call cfwure (nput, nFc*2,    st % dv % fch_52, 'fch_52')
      call cfwure (nput, nCv*ns,   st % dv % kinrgy, 'kinrgy')
      call cfwure (nput, nFc*2*ns, st % dv % fhm, 'fhm')
      call cfwure (nput, nFc*2,    st % dv % fkt, 'fkt')

*   ..save time
      time(1)=tim
      call cfwure (nput, 1, time, 'time')

*   ..save parallel current
      call cfwure (nput, nFc*2,    st % dv % fch_p, 'fch_p')

*   ..saving various fluxes and velocities needed by advanced SOLPS-ITER features
      call cfwure (nput, nFc*2*ns, st % dv % fna_mdf, 'fna_mdf')
      call cfwure (nput, nFc*2,    st % dv % fhe_mdf, 'fhe_mdf')
      call cfwure (nput, nFc*2,    st % dv % fhi_mdf, 'fhi_mdf')
      call cfwure (nput, nFc*2*ns, st % dv % fna_fcor, 'fna_fcor')
      call cfwure (nput, nFc*2*ns, st % dv % fna_nodrift, 'fna_nodrift')
      call cfwure (nput, nFc*2*ns, st % dv % fna_he, 'fna_he')
      call cfwure (nput, nFc*2*ns, st % dv % fnaPSch, 'fnaPSch')
      call cfwure (nput, nFc*2,    st % dv % fhePSch, 'fhePSch')
      call cfwure (nput, nFc*2,    st % dv % fhiPSch, 'fhiPSch')
      call cfwure (nput, nFc*2*ns, st % dv % fna_eir, 'fna_eir')
      call cfwure (nput, nFc*2,    st % dv % fne_eir, 'fne_eir')
      call cfwure (nput, nFc*2,    st % dv % fhe_eir, 'fhe_eir')
      call cfwure (nput, nFc*2,    st % dv % fhi_eir, 'fhi_eir')
      call cfwure (nput, nFc*2*ns, st % dv % fna_32, 'fna_32')
      call cfwure (nput, nFc*2*ns, st % dv % fna_52, 'fna_52')
      call cfwure (nput, nFc*2,    st % dv % fni_32, 'fni_32')
      call cfwure (nput, nFc*2,    st % dv % fni_52, 'fni_52')
      call cfwure (nput, nFc*2,    st % dv % fne_32, 'fne_32')
      call cfwure (nput, nFc*2,    st % dv % fne_52, 'fne_52')
      call cfwure (nput, nFc*2,    st % dv % fchdia, 'fchdia')
      call cfwure (nput, nFc*2,    st % dv % fchin, 'fchin')
      call cfwure (nput, nFc*2,    st % dv % fchvispar, 'fchvispar')
      call cfwure (nput, nFc*2,    st % dv % fchvisper, 'fchvisper')
      call cfwure (nput, nFc*2,    st % dv % fchvisq, 'fchvisq')
      call cfwure (nput, nFc*2,    st % dv % fchinert, 'fchinert')
      call cfwure (nput, nFc*2*ns, st % dv % vaecrb, 'vaecrb')
      call cfwure (nput, nFc*2*ns, st % dv % vadia, 'vadia')
      call cfwure (nput, nFc*2*ns, st % dv % wadia, 'wadia')
      call cfwure (nput, nFc*2,    st % dv % veecrb, 'veecrb')
      call cfwure (nput, nFc*2,    st % dv % vedia, 'vedia')
      call cfwure (nput, nFc*2,    st % dv % floe_noc, 'floe_noc')
      call cfwure (nput, nFc*2,    st % dv % floi_noc, 'floi_noc')

* ..return
      call subend ()
      return

      end subroutine write_b2fstate_diff

