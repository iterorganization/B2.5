      subroutine b2xddr (ninp, nout, nx, ny, ns)
      use b2mod_types
      use b2mod_geo
      use b2mod_work
      use b2mod_plasma
      use b2mod_rates
      use b2mod_residuals
      use b2mod_sources
      use b2mod_transport
      use b2mod_indirect
      use b2mod_constants
      use b2mod_geo_corner
      use b2mod_b2cmgs
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_b2cmpt
      use b2mod_b2cmwg
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer ninp(0:4), nout(0:2), nx, ny, ns
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
*   ..local variables
      integer is, idum(0:9)
      real (kind=R8) :: za(0:ns-1)
      character lblgm*120, lblcp*120, lblmn*120
*   ..procedures
      external xertst
      external cfruch, cfruin
      external b2rugm, b2rucp, b2ruzd, b2rups, b2xvps, process

*   ..subprogram start-up calls
      call subini ('b2xddr')
*   ..test ninp, nout
      call xertst (1.le.ninp(0).and.1.le.ninp(1).and.1.le.ninp(2).and.
     &             1.le.ninp(3).and.1.le.nout(0).and.1.le.nout(1).and.
     &             1.le.nout(2),
     &  'faulty argument ninp, nout')
*   ..test dimensions
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (ns.le.nsdecl, 'faulty parameter nsdecl')

* ..computation
*   ..read geometry
      call cfruch (ninp(1), 120, lblgm, 'label')
      call cfruin (ninp(1), 1, idum, 'isymm')
      call b2rugm (ninp(1), nx, ny, crx, cry, fpsi, ffbz,
     &  bb, vol, hx, hy, qz, qc, qs, gs, pbs,
     &  wbbl, wbbr, wbbv, wbbc,
     &  cell_width, cell_height, gmap)
*   ..read physics common
      call cfruch (ninp(2), 120, lblcp, 'label')
      call b2rucp (ninp(2), b2fparam_version, nx, ns)
*   ..read plasma state file
      call cfruch (ninp(3), 120, lblmn, 'label')
*   ..read species data
      call b2ruzd (ninp(3), b2fstate_version,
     &  ns, zamin, zamax, zn, am, .true.)
      do is = 0, ns-1
        za(is) = (zamin(is) + zamax(is))/2.0_R8
      end do
*   ..read plasma state
      call b2rups (ninp(3), nx, ny, ns, ne, na, ua, uadia, te,          !xpb
     &  ti, po, fna, fhe, fhi, fch, fch_32, fch_52, fch_p, kinrgy,
     &  fna_mdf, fhe_mdf, fhi_mdf, fna_fcor, fna_nodrift, fna_he,
     &  fnaPSch, fhePSch, fhiPSch, fna_eir, fne_eir, fhe_eir, fhi_eir,
     &  fna_32, fna_52, fni_32, fni_52, fne_32, fne_52,
     &  fchdia, fchin, fchvispar, fchvisper, fchvisq, fchinert,
     &  vaecrb, vadia, wadia, veecrb, vedia, floe_noc, floi_noc,
     &  b2fstate_version)
*   ..test plasma state and fluxes
      call b2xvps (nx, ny, ns,
     &  ne, na, ua, te, ti, po, fna, fhe, fhi, fch)
*   ..process
      call process(nx,ny,ns,
     1  crx,cry,bb,vol,hx,hy,qz,gs,
     2  za,zn,am,na,ua,uadia,te,ti,po,fna,fhe,fhi,fch,fchdia,
     3  resco,reshe,reshi,resmo,resmt,respo,
     4  nncut,leftcut,rightcut,topcut,bottomcut)

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2xddr

      end subroutine b2xddr
