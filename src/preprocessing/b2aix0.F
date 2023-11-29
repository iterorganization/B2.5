*-----------------------------------------------------------------------
*.specification

      subroutine b2aix0 (ninp, nCv, nFc, ns, mpg)
      use b2mod_types
      use b2mod_subsys
      use b2mod_version
      use b2us_map
      implicit none
*   ..input arguments (unchanged on exit)
      integer ninp(0:2)
*   ..output arguments (unspecified on entry)
      integer nCv, nFc, ns
      type(mapping), intent(inout) :: mpg

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2AIX0 reads the mesh dimensions and the number of species from a
*     formatted file.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine reads from unit ninp(0) the number of species, ns,
*     and from ninp(1) the mesh dimensions. These units
*     must be connected to a formatted file. The input files must be
*     positioned properly: the first record must begin with the text
*     '*dimens ' or '*DIMENS ', and the next record must contain the
*     input values in free format. Only those two records are read
*     from each of the two input files.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     ninp - (0:2) integer array.
*     ninp specifies the input unit numbers.
*
*     nCv, nFc - integer.
*     nCv and nFc specify the number of control volumes and cell faces
*     in the unstructured mesh.
*     It will hold that 0.le.nCv and 0.le.nFc.
*
*     ns - integer, output.
*     ns specifies the number of atomic species in the calculation.
*     The species are indexed by (0:ns-1).
*     It will hold that 1.le.ns.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer idum(0:6)
      character id*8
*   ..procedures
      logical streql
      external xertst, xerrab, streql

*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2aix0')
*   ..test ninp
      call xertst (1.le.ninp(0).and.1.le.ninp(1),
     &  'faulty argument ninp')
*   ..read header line for ns
      read (ninp(0),'(a8)',err=93) id
      call xertst (streql(id,'*dimens'),
     &  'faulty input id at *dimens for ns; id='//id)
*   ..read ns
      read (ninp(0),*,err=94) ns
      call xertst (1.le.ns, 'faulty input ns')
*   ..obtain mapping/geometry dimensions from file
      call cfverr (ninp(1), b2fgmtry_version)
      if (b2fgmtry_version .lt. '03.002.000') then
        call xerrab('Using old input format. Run b2us first.')
      endif
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
*   ..test basic dimensions
      nCv = mpg%nCv
      nFc = mpg%nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty input nCv, nFc')

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.error conditions

   93 call xerrab ('error trying to read heading for ns')
   94 call xerrab ('error trying to read ns')

*-----------------------------------------------------------------------
*.end b2aix0

      end subroutine b2aix0
