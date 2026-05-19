*-----------------------------------------------------------------------
*.specification

      program b2ch
      use b2mod_types
      use b2mod_version
      use b2mod_xerset
      use b2us_data
      use b2us_geo
      use b2us_map
      implicit none

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2CH (main program) prepares input files for the ECH module in B2.
*     This program calls the driver routine b2chdr for the actual
*     computation.
*
*
*  2. specification
*
*     Main program.
*
*
*  3. description (see also routine b2cdca)
*
*     The complete program b2ch generates the chords which mimic the
*     path of ECH rays, taking into account their bouncing on vessel
*     walls. It also provides the location of the resonance layer and
*     the control volumes crossed by the rays, in which the ECH power
*     will be deposited.
*     This program unit computes and writes to external files the chords
*     coordinates and their respective labels, the location of the ECH
*     resonance layer and a list of control volumes.
*
*     This program unit opens and closes the input/output units, and
*     may perform some other system-dependent operations. It calls
*     routine b2chdr to control the actual computation.
*
*     The input units are:
*     ninp(0): formatted; provides chords parameters.
*     ninp(1): un*formatted; provides the geometry.
*
*     The output units are:
*     nout(0): formatted; provides printed output.
*     nout(1): formatted; provides the chords and their labels.
*     nout(2): un*formatted; provides the location of the resonance
*              layer, the list of control volumes and the relative
*              fraction of deposited power. !!TODO
*
*     (See routine b2cdca for the meaning of 'un*formatted'.)
*
*
*  4. references
*
*     ...
*
*
*  5. parameters (see also routine b2cdcv)
*
*     None.
*
*
*  6. error indicators
*
*     In case an error condition is detected, a call is made to the
*     routine xerrab. This causes an error message to be printed,
*     after which the program halts.
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks

*   ..local variables
      integer ninp(0:2), nout(0:2)
      integer ns
      integer io, idum(0:9)
      logical file_open
*   ..procedures
      character*32 get_B25_hash
      external get_B25_hash
      external prgini, prgend, xertst, cfopen, cfruin
      external b2chdr
*   ..initialise input/output units
      data ninp/50,51,52/
      data nout/60,61,62/

*-----------------------------------------------------------------------
*.documentation-internal
*
*     The following common blocks have their outermost declaration in
*     this routine; they need not be preserved between calls.
*
*     Description of some local variables:
*
*     ninp - (0:1) integer array.
*     ninp specifies the input unit numbers.
*
*     nout - (0:2) integer array.
*     nout specifies the output unit numbers.
*
*-----------------------------------------------------------------------
*.computation

*   ..program start-up calls
      call prgini ('b2ch')
*   ..open input files
      call cfopen (ninp(0),'b2ch.dat','old','formatted')
      call cfopen (ninp(1),'b2fgmtry','old','un*formatted')
      call cfopen (ninp(2),'structure.dat','old','formatted')
*   ..mark output unit for error messages
      call xerset (0)
*   ..obtain version numbers
      call cfverr (ninp(1), b2fgmtry_version)
      write(6,*) 'Running with VERSION ',newversion,' ', get_B25_hash()
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
*   ..obtain ns from parameters file
*   ..call principal driver routine
      call alloc_mapping (mpg)
      call alloc_geometry (geo, mpg%nCv, mpg%nFc, mpg%nVx, mpg%nFs,
     .  mpg%nFt, mpg%nVmxCv)
*   ..open output files
      call cfopen (nout(0),'b2ch.prt','new','formatted')
      call cfopen (nout(1),'rays.chords','new','formatted')
      call cfopen (nout(2),'b2chords','new','un*formatted')
*   ..compute, display and save chords
      call b2chdr (ninp, nout, ns)
*   ..close files
      do io = 0, size(ninp)-1
        inquire(unit=ninp(io),opened=file_open)
        if (file_open) close (ninp(io))
      enddo
      do io = 0, size(nout)-1
        inquire(unit=nout(io),opened=file_open)
        if (file_open) close(nout(io))
      enddo

* ..end of computation
      call dealloc_geometry (geo)
      call dealloc_mapping (mpg)
      call prgend ()
      stop 'b2ch'

*-----------------------------------------------------------------------
*.end b2ch

      end program b2ch
