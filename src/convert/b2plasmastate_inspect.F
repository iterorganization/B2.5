      program b2plasmastate_inspect
      use b2mod_types
      use b2mod_version
      use b2mod_time
      use b2mod_geo
      use b2mod_plasma
      use b2mod_rates
      use b2mod_residuals
      use b2mod_sources
      use b2mod_sputter
      use b2mod_transport
      use b2mod_anomalous_transport
      use b2mod_work
      use b2mod_indirect
#ifdef NAGF95
      use f90_unix           ! IGNORE
#endif
      use b2mod_xerset
      implicit none
#include "DIMENSIONS.F"

*-----------------------------------------------------------------------
*     b2plasmastate_inspect reads in a formatted plasmastate file, and
*     displays the value of the plasma time stamp "tim"
*
*     Run as: b2plasmastate_inspect.exe plasmastate.xxxx
*-----------------------------------------------------------------------
      
#ifdef USE_PXFGETENV
      integer lenval, ierror
#else
#ifdef NAGFOR
      integer lenval, ierror
#endif
#endif
#ifndef NAGF95
#ifndef NAGFOR
      integer iargc
#endif
#endif
      integer :: ninp(0:5), nx, ny, ns, ns0, idum(0:9)
      character :: arg*256
*   ..procedures
      external xertst, cfopen, cfruin
      external b2rucp, b2rups
*   ..initialise input/output unit numbers
      data ninp/50,51,52,53,54,55/

*-----------------------------------------------------------------------

      call prgini('b2plasmastate_inspect')
! Handle command args
#ifndef NAGFOR
      call xertst (iargc().ge.1,
     &  'program should be supplied with a state number')
#ifdef USE_PXFGETENV
      call pxfgetarg(1,arg,lenval,ierror)
#else
      call getarg(1,arg)
#endif
#else
      call xertst (command_argument_count().ge.1,
     &  'program should be supplied with a state number')
      call get_command_argument(1,arg,lenval,ierror)
#endif

! Open files
      call cfopen (ninp(1),'b2fgmtry','old','un*formatted')
      call cfopen (ninp(2),'b2fpardf','old','un*formatted')
      call cfopen (ninp(3),'b2fstati','old','un*formatted')
      call cfopen (ninp(4), trim(arg), 'old','unformatted')

! Follow b2co logic      
*   ..mark output unit for error messages
      call xerset (0)
*   ..obtain version numbers                                           !xpb
      call cfverr (ninp(1), b2fgmtry_version)                          !xpb
      call cfverr (ninp(2), b2fparam_version)                          !xpb
      call cfverr (ninp(3), b2fstati_version)                          !xpb
*   ..obtain nx, ny from geometry file
      call cfruin (ninp(1), 2, idum, 'nx,ny')
      nx = idum(0)
      ny = idum(1)
      call xertst (0.le.nx.and.0.le.ny,
     &  'faulty input nx, ny from geometry file')
*   ..obtain ns from parameters file
      call cfruin (ninp(2), 1, idum, 'ns')
      ns = idum(0)
      call xertst (1.le.ns, 'faulty input ns from parameters file')
*   ..obtain nx, ny, ns0 from plasma state file
      call cfruin (ninp(3), 3, idum, 'nx,ny,ns')
      ns0 = idum(2)
      call xertst (idum(0).eq.nx.and.idum(1).eq.ny,
     &  'faulty input nx, ny from plasma state file')
      call xertst (1.le.ns0.and.ns0.le.ns,
     &  'faulty input ns from plasma state file')
*   ..call principal driver routine
      call alloc_b2mod_geo(nx,ny)
      call alloc_b2mod_plasma(nx,ny,ns)
      call alloc_b2mod_rates(nx,ny,ns)
      call alloc_b2mod_residuals(nx,ny,ns)
      call alloc_b2mod_sources(nx,ny,ns)
      call alloc_b2mod_sputter(nx,ny,ns)
      call alloc_b2mod_transport(nx,ny,ns)
      call alloc_b2mod_anomalous_transport(nx,ny,ns)
      call alloc_b2mod_work(nx,ny,ns)
      call alloc_b2mod_indirect(nx,ny,nncutmax)                         !xpb

*     ..read initial state
      call cfverr (ninp(4), plasmastate_version)                        !xpb
      call b2rups (ninp(4), nx, ny, ns, ne, na, ua, uadia, te, ti,      !xpb
     &  po, fna, fhe, fhi, fch, fch_32, fch_52, fch_p, kinrgy,
     &  fna_mdf, fhe_mdf, fhi_mdf, fna_fcor, fna_nodrift, fna_he,
     &  fnaPSch, fhePSch, fhiPSch, fna_eir, fne_eir, fhe_eir, fhi_eir,
     &  fna_32, fna_52, fni_32, fni_52, fne_32, fne_52,
     &  fchdia, fchin, fchvispar, fchvisper, fchvisq, fchinert,
     &  vaecrb, vadia, wadia, veecrb, vedia, floe_noc, floi_noc,
     1  plasmastate_version)

! Display tim      
      write(0,'(3a,1p,1e15.8,a)') 'The plasmastate file ',trim(arg),
     & ' corresponds to plasma time ',tim,' seconds.'

* ..end of computation
      call dealloc_b2mod_geo
      call dealloc_b2mod_work
      call dealloc_b2mod_rates
      call dealloc_b2mod_plasma
      call dealloc_b2mod_sources
      call dealloc_b2mod_sputter
      call dealloc_b2mod_indirect
      call dealloc_b2mod_residuals
      call dealloc_b2mod_transport
      call dealloc_anomalous_transport
      call prgend ()

      stop
      end program b2plasmastate_inspect

