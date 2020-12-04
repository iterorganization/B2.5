*-----------------------------------------------------------------------
*.specification

      subroutine b2aix0 (ninp, nx, ny, ns)
      use b2mod_types
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer ninp(0:1)
*   ..output arguments (unspecified on entry)
      integer nx, ny, ns

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
*     This routine reads from unit ninp(0) the mesh dimensions, nx and
*     ny, and from unit ninp(1) the number of species, ns. These units
*     must be connected to a formatted file. The input files must be
*     positioned properly: the first record must begin with the text
*     '*dimens ' or '*DIMENS ', and the next record must contain the
*     input values in free format. Only those two records are read
*     from each of the two input files.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     ninp - (0:1) integer array.
*     ninp specifies the input unit numbers.
*
*     nx, ny - integer, output.
*     nx and ny specify the number of interior cells along the first
*     and the second coordinate, respectively. The total number of
*     cells is (nx+2)*(ny+2); they are indexed by (-1:nx,-1:ny).
*     It will hold that 0.le.nx and 0.le.ny.
*
*     ns - integer, output.
*     ns specifies the number of atomic species in the calculation.
*     The species are indexed by (0:ns-1).
*     It will hold that 1.le.ns.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
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
*   ..read header line for nx, ny
      read (ninp(0),'(a8)',err=91) id
      call xertst (streql(id,'*dimens'),
     &  'faulty input id at *dimens for nx, ny; id='//id)
*   ..read nx, ny
      read (ninp(0),*,err=92) nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty input nx, ny')
*   ..read header line for ns
      read (ninp(1),'(a8)',err=93) id
      call xertst (streql(id,'*dimens'),
     &  'faulty input id at *dimens for ns; id='//id)
*   ..read ns
      read (ninp(1),*,err=94) ns
      call xertst (1.le.ns, 'faulty input ns')

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.error conditions

   91 call xerrab ('error trying to read heading for nx, ny')
   92 call xerrab ('error trying to read nx, ny')
   93 call xerrab ('error trying to read heading for ns')
   94 call xerrab ('error trying to read ns')

*-----------------------------------------------------------------------
*.end b2aix0

      end subroutine b2aix0
