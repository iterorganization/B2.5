*-----------------------------------------------------------------------
*.specification

      program b2fu
      use b2mod_types
      use b2mod_version
      use b2mod_b2cmpa
      use b2mod_xerset
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2us_io
      implicit none

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2FU reads in a formatted b2fplasmf file, and writes out a
*     an unformatted b2fplasma file.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks

*   ..main data structures
      type(mapping), save  :: mpg
      type(geometry), save :: geo
      type(B2State), save  :: state
*   ..local variables
      integer ninp(0:5), nout(0:7), nCv, nFc, ns, idum(0:9), io
      logical file_open
      character(120)    :: cdum
*   ..procedures
      external prgini, prgend, xertst, cfopen, cfruin
*   ..initialise input/output unit numbers
      data ninp/50,51,52,53,54,55/
      data nout/60,61,62,63,64,65,66,67/

*-----------------------------------------------------------------------
*.computation

*   ..program start-up calls
      call prgini ('b2fu')
*   ..open files
      call cfopen (ninp(1),'b2fgmtry','old','un*formatted')
      call cfopen (ninp(2),'b2fparam','old','un*formatted')
      call cfopen (ninp(3),'b2fstate','old','un*formatted')
      call cfopen (ninp(4),'b2fplasmf','old','formatted')
      call cfopen (nout(0),'b2fu.prt','new','formatted')
      call cfopen (nout(7),'b2fplasma','new','unformatted')
*   ..mark output unit for error messages
      call xerset (0)
*   ..obtain version numbers                                             !xpb
      call cfverr (ninp(1), b2fgmtry_version)                            !xpb
      call cfverr (ninp(2), b2fparam_version)                            !xpb
      call cfverr (ninp(3), b2fstate_version)                            !xpb
      call cfverr (ninp(4), b2fplasma_version)
*   ..obtain nCv, nFc, ns from b2fstate file
      call cfruin (ninp(3), 3, idum, 'nCv,nFc,ns')
      nCv = idum(0)
      nFc = idum(1)
      ns  = idum(2)
      call xertst (0.lt.nCv.and.0.lt.nFc.and.1.le.ns,
     &  'faulty input nCv, nFc, ns from plasma state file')
*   ..obtain ns from b2fparam file
      call cfruin (ninp(2), 1, idum, 'ns')
      call xertst (idum(0).eq.ns,
     &  'faulty input ns from physics parameters file')
*   ..obtain mapping/geometry dimensions from b2fgmtry file
      call read_b2fgmtry_header(ninp(1), b2fgmtry_version, 
     &  mpg, cdum, idum(0))
      call xertst (mpg%nCv.eq.nCv.and.mpg%nFc.eq.nFc,
     &  'faulty input nCv, nFc from geometry file')
*   ..obtain parameters from b2fplasmf file
      call cfruin (ninp(4),3,idum,'nCv,nFc,ns')
      call xertst (idum(0).eq.nCv.and.idum(1).eq.nFc.and.idum(2).eq.ns,
     &  'faulty input nCv, nFc, ns from b2fplasmf file')
*   ..call principal driver routine
      call alloc_mapping (mpg)
      call alloc_geometry (geo, mpg%nCv, mpg%nFc, mpg%nVx, mpg%nFs,
     .  mpg%nFt, mpg%nVmxCv)
      call createB2State (nCv, nFc, ns, state)
*   ..read plasma state
      call read_b2fplasma(ninp(4), nCv, nFc, ns, state)
*   ..write plasma state
      newversion=b2fplasma_version
      call cfverw(nout(7), newversion)
      idum(0) = nCv
      idum(1) = nFc
      idum(2) = ns
      call cfwuin (nout(7), 3, idum, 'nCv,nFc,ns')
      call write_b2fplasma (nout(7), nCv, nFc, ns, state)
*   ..close files
      do io=0,size(ninp)-1
         inquire(unit=ninp(io),opened=file_open)
         if (file_open) close(ninp(io))
      enddo
      do io=0,size(nout)-1
         inquire(unit=nout(io),opened=file_open)
         if (file_open) close(nout(io))
      enddo
* ..end of computation
      call destroyB2State (state)
      call dealloc_geometry (geo)
      call dealloc_mapping (mpg)
      call prgend ()
      stop 'b2fu'

*-----------------------------------------------------------------------
*.end b2fu

      end program b2fu
