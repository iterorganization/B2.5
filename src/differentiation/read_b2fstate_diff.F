      subroutine read_b2fstate_diff (nput, nCv, nFc, ns, st)
      use b2mod_time
      use b2us_plasma_diff
      implicit none
      integer, intent(in) :: nput, nCv, nFc, ns
      type(B2STATE_DIFF), intent(inout) :: st
      integer :: ix, iy, ic, ifc, nfx, idum(3)
      real (kind=R8) :: time(1)

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('read_b2fstate_diff')
*   ..test nput
      call xertst (1.le.nput, 'faulty argument nput')
*   ..test ns
      call xertst (1.le.ns, 'faulty argument ns')
!csc      call createB2State (nCv, nFc, ns, st)
     
* ..read state and fluxes

*   ..read plasma state
      call cfrure (nput, nCv*ns,   st % pl % na, 'na')
      call cfrure (nput, nCv,      st % dv % ne, 'ne')
      call cfrure (nput, nCv*ns,   st % pl % ua, 'ua')
      call cfrure (nput, nFc*2*ns, st % dv % uadia, 'uadia') 
      call cfrure (nput, nCv,      st % pl % te, 'te')
      call cfrure (nput, nCv,      st % pl % ti, 'ti')
      call cfrure (nput, nCv,      st % pl % tn, 'tn')
      call cfrure (nput, nCv,      st % pl % po, 'po')
      call cfrure (nput, nCv,      st % pl % kt, 'kt')

*   ..read fluxes
      call cfrure (nput, nFc*2*ns, st % dv % fna, 'fna')
      call cfrure (nput, nFc*2,    st % dv % fhe, 'fhe')
      call cfrure (nput, nFc*2,    st % dv % fhi, 'fhi')
      call cfrure (nput, nFc*2,    st % dv % fhn, 'fhn')
      call cfrure (nput, nFc*2,    st % dv % fch, 'fch')
      call cfrure (nput, nFc*2,    st % dv % fch_32, 'fch_32')
      call cfrure (nput, nFc*2,    st % dv % fch_52, 'fch_52')
      call cfrure (nput, nCv*ns,   st % dv % kinrgy, 'kinrgy')
      call cfrure (nput, nFc*2*ns, st % dv % fhm, 'fhm')
      call cfrure (nput, nFc*2,    st % dv % fkt, 'fkt')

*   ..read time
      call cfrure (nput, 1, time, 'time')
      tim = time(1)

*   ..read parallel current
      call cfrure (nput, nFc*2,    st % dv % fch_p, 'fch_p')

*   ..read various fluxes and velocities needed by advanced SOLPS-ITER features
      call cfrure (nput, nFc*2*ns, st % dv % fna_mdf, 'fna_mdf')
      call cfrure (nput, nFc*2,    st % dv % fhe_mdf, 'fhe_mdf')
      call cfrure (nput, nFc*2,    st % dv % fhi_mdf, 'fhi_mdf')
      call cfrure (nput, nFc*2*ns, st % dv % fna_fcor, 'fna_fcor')
      call cfrure (nput, nFc*2*ns, st % dv % fna_nodrift, 'fna_nodrift')
      call cfrure (nput, nFc*2*ns, st % dv % fna_he, 'fna_he')
      call cfrure (nput, nFc*2*ns, st % dv % fnaPSch, 'fnaPSch')
      call cfrure (nput, nFc*2,    st % dv % fhePSch, 'fhePSch')
      call cfrure (nput, nFc*2,    st % dv % fhiPSch, 'fhiPSch')
      call cfrure (nput, nFc*2*ns, st % dv % fna_eir, 'fna_eir')
      call cfrure (nput, nFc*2,    st % dv % fne_eir, 'fne_eir')
      call cfrure (nput, nFc*2,    st % dv % fhe_eir, 'fhe_eir')
      call cfrure (nput, nFc*2,    st % dv % fhi_eir, 'fhi_eir')
      call cfrure (nput, nFc*2*ns, st % dv % fna_32, 'fna_32')
      call cfrure (nput, nFc*2*ns, st % dv % fna_52, 'fna_52')
      call cfrure (nput, nFc*2,    st % dv % fni_32, 'fni_32')
      call cfrure (nput, nFc*2,    st % dv % fni_52, 'fni_52')
      call cfrure (nput, nFc*2,    st % dv % fne_32, 'fne_32')
      call cfrure (nput, nFc*2,    st % dv % fne_52, 'fne_52')
      call cfrure (nput, nFc*2,    st % dv % fchdia, 'fchdia')
      call cfrure (nput, nFc*2,    st % dv % fchin, 'fchin')
      call cfrure (nput, nFc*2,    st % dv % fchvispar, 'fchvispar')
      call cfrure (nput, nFc*2,    st % dv % fchvisper, 'fchvisper')
      call cfrure (nput, nFc*2,    st % dv % fchvisq, 'fchvisq')
      call cfrure (nput, nFc*2,    st % dv % fchinert, 'fchinert')
      call cfrure (nput, nFc*2*ns, st % dv % vaecrb, 'vaecrb')
      call cfrure (nput, nFc*2*ns, st % dv % vadia, 'vadia')
      call cfrure (nput, nFc*2*ns, st % dv % wadia, 'wadia')
      call cfrure (nput, nFc*2,    st % dv % veecrb, 'veecrb')
      call cfrure (nput, nFc*2,    st % dv % vedia, 'vedia')
      call cfrure (nput, nFc*2,    st % dv % floe_noc, 'floe_noc')
      call cfrure (nput, nFc*2,    st % dv % floi_noc, 'floi_noc')

* ..return
      call subend ()
      return

      end subroutine read_b2fstate_diff

