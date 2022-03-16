*-----------------------------------------------------------------------
*.specification

      subroutine b2ahx0 (nrid, ns)
      use b2mod_types
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nrid
*   ..output arguments (unspecified on entry)
      integer ns

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2AHX0 reads the number of species from a formatted file.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine reads from unit nrid the number of species, ns.
*     Unit nrid must be connected to a formatted file. The input file
*     must be positioned properly: the first record must begin with the
*     text '*dimens ' or '*DIMENS ', and the next record must contain
*     the value of ns in free format. Only those two records are read.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nrid - integer, input.
*     nrid specifies a unit number for formatted input.
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
      call subini ('b2ahx0')
*   ..test nrid
      call xertst (1.le.nrid, 'faulty argument nrid')
*   ..read header line
      read (nrid,'(a8)',err=91) id
      call xertst (streql(id,'*dimens'),
     &  'faulty input id at *dimens, id='//id)
*   ..read ns
      read (nrid,*,err=92) ns
      call xertst (1.le.ns, 'faulty input ns')

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.error conditions

   91 call xerrab ('error trying to read heading *dimens')
   92 call xerrab ('error trying to read after '//id)

*-----------------------------------------------------------------------
*.end b2ahx0

      end
