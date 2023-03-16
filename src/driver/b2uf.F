*-----------------------------------------------------------------------
*.specification

      program b2uf
      use b2mod_types
      use b2mod_version
      use b2mod_b2cmpa
      use b2mod_b2cmfs
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2us_prep
      use b2us_io
      use b2mod_xerset
      implicit none

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2UF reads in an unformatted b2fplasma file, and writes out a
*     a formatted b2fplasma (b2fplasmf) ile.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks

*   ..main data structures
      type(mapping), save    :: mpg
      type(geometry), save   :: geo
      type(B2State), save    :: state
!      type(B2StateExt), save :: state_ext
!      type(switches), save   :: switch

*   ..local variables
      integer ninp(0:5), nout(0:7), nCv, nFc, ns, idum(0:9), io
      logical file_open
      character lblgm*120, lblps*120
*   ..procedures
      external prgini, prgend, xertst, cfopen, cfruin
*   ..initialise input/output unit numbers
      data ninp/50,51,52,53,54,55/
      data nout/60,61,62,63,64,65,66,67/

*-----------------------------------------------------------------------
*.computation

*   ..program start-up calls
      call prgini ('b2uf')
*   ..open files
      call cfopen (ninp(1),'b2fgmtry','old','un*formatted')
!      call cfopen (ninp(2),'b2fparam','old','un*formatted')
      call cfopen (ninp(3),'b2fstate','old','un*formatted')
      call cfopen (ninp(4),'b2fplasma','old','unformatted')
      call cfopen (nout(0),'b2uf.prt','new','formatted')
      call cfopen (nout(7),'b2fplasmf','new','un*formatted')
*   ..mark output unit for error messages
      call xerset (0)
*   ..obtain version numbers
      call cfverr (ninp(1), b2fgmtry_version)
!      call cfverr (ninp(2), b2fparam_version)
      call cfverr (ninp(3), b2fstate_version)
      call cfverr (ninp(4), b2fplasma_version)
*   ..obtain nCv, nFc, ns from b2fstate file
      call cfruin (ninp(3), 3, idum, 'nCv,nFc,ns')
      nCv = idum(0)
      nFc = idum(1)
      ns  = idum(2)
      call xertst (0.le.nCv.and.0.le.nFc.and.1.le.ns,
     &  'faulty input nCv, nFc, ns from b2fstate file')
*   ..obtain ns from parameters file
!      call cfruin (ninp(2), 1, idum(0), 'ns')
!      call xertst (idum(0).eq.ns,
!     &  'faulty input ns from physics parameters file')
*   ..obtain mapping/geometry dimensions from b2fgmtry file
      call cfruin (ninp(1),7,idum,'nCi,nCg,nCv,nFc,nVx,nFs,nFt')
      mpg%nCi = idum(0)
      mpg%nCg = idum(1)
      mpg%nCv = idum(2)
      mpg%nFc = idum(3)
      mpg%nVx = idum(4)
      mpg%nFs = idum(5)
      mpg%nFt = idum(6)
      call cfruin (ninp(1),5,idum,'nCmxVx,nCmxFc,nVmxCv,nVmxFc,nCmxNv')
      mpg%nCmxVx = idum(0)
      mpg%nCmxFc = idum(1)
      mpg%nVmxCv = idum(2)
      mpg%nVmxFc = idum(3)
      mpg%nCmxNv = idum(4)
      call cfruin (ninp(1),1,idum,'isClassicalGrid')
      mpg%isClassicalGrid = idum(0)
      call cfruin (ninp(1),3,idum,'nx,ny,nncut')
      mpg%nx    = idum(0)
      mpg%ny    = idum(1)
      mpg%nncut = idum(2)
      call xertst (mpg%nCv.eq.nCv.and.mpg%nFc.eq.nFc,
     &  'faulty input nCv, nFc from geometry file')
*   ..obtain parameters from b2fplasma file
      call cfruin (ninp(4), 3, idum, 'nCv,nFc,ns')
      call xertst (idum(0).eq.nCv.and.idum(1).eq.nFc.and.idum(2).eq.ns,
     &  'faulty input nCv, nFc, ns from b2fplasma file')
*   ..call principal driver routine
      call alloc_mapping (mpg)
      call alloc_geometry (geo, mpg%nCv, mpg%nFc, mpg%nVx, mpg%nFs,
     .  mpg%nFt, mpg%nVmxCv)
      call createB2State (nCv, nFc, ns, state)
!      call createB2StateExt ( mpg%nCv, mpg%nFc, ns, state_ext )
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
*   ..for underlying structured grids, convert to structured output
      if (mpg%isClassicalGrid.eq.1) then
*    ..read header for geometry
       call cfruch (ninp(1), 120, lblgm, 'label')
#ifdef JET
       labgeo=lblgm
#endif
       call cfruin (ninp(1), 1, idum, 'isymm')
       isymm = idum(0)
*      (isymm is now stored in the b2cmfs common)
*    ..read header for initial state
*      (zamin, zamax, zn, am from ninp(3) will be checked against the data
*      already read in)
       call cfruch (ninp(3), 120, lblps, 'label')
       call b2ruzd (ninp(3), b2fstate_version, ns, 
     &  zamin, zamax, zn, am, .false.)
*    ..read full b2fstate and b2fgmtry files
       call read_b2fgmtry (ninp(1), mpg, geo)
       call read_b2fstate (ninp(3), mpg%nCv, mpg%nFc, ns, state)
*    ..convert to structured format and write b2fstate_st, b2fplasmf_st
       call to_struct_plasma (ns,mpg,state)
      endif
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
!      call destroyB2StateExt (state_ext)
      call dealloc_geometry (geo)
      call dealloc_mapping (mpg)
      call prgend ()
      stop 'b2uf'

*-----------------------------------------------------------------------
*.end b2uf

      end program b2uf
