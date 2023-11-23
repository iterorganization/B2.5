      subroutine b2xddr (ninp, nout, ns, mpg, geo, state)
      use b2mod_types
      use b2mod_b2cmgs
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_b2cmpt
      use b2mod_b2cmwg
      use b2mod_subsys
      use b2mod_version
      use b2us_io
      use b2us_geo
      use b2us_map
      use b2us_plasma
      implicit none
*   ..input arguments (unchanged on exit)
      integer ninp(0:4), nout(0:2), ns
*   ..main data structures
      type(mapping)  :: mpg
      type(geometry) :: geo
      type(B2State)  :: state
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
      call xertst (0.le.mpg%nCv.and.0.le.mpg%nFc,
     &  'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (ns.le.nsdecl, 'faulty parameter nsdecl')

* ..computation
*   ..read geometry
      call cfruch (ninp(1), 120, lblgm, 'label')
      call cfruin (ninp(1), 1, idum, 'isymm')
      call read_b2fgmtry (ninp(1), mpg, geo)
*   ..read physics common
      call cfruch (ninp(2), 120, lblcp, 'label')
      call b2rucp (ninp(2), b2fparam_version, ns)
*   ..read plasma state file
      call cfruch (ninp(3), 120, lblmn, 'label')
*   ..read species data
      call b2ruzd (ninp(3), b2fstate_version,
     &  ns, zamin, zamax, zn, am, .true.)
      do is = 0, ns-1
        za(is) = (zamin(is) + zamax(is))/2.0_R8
      end do
*   ..read plasma state
      call read_b2fstate (ninp(3), mpg%nCv, mpg%nFc, ns, state)
*   ..test plasma state and fluxes
      call b2xvps (mpg%nCv, mpg%nFc, ns, state%pl, state%dv)
*   ..process
#ifdef WG_TODO
      call process(nx,ny,ns,
     1  crx,cry,bb,vol,hx,hy,qz,gs,
     2  za,zn,am,na,ua,uadia,te,ti,po,fna,fhe,fhi,fch,fchdia,
     3  resco,reshe,reshi,resmo,resmt,respo,
     4  nncut,leftcut,rightcut,topcut,bottomcut)
#endif

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2xddr

      end subroutine b2xddr
