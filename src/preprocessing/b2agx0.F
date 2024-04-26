*-----------------------------------------------------------------------
*.specification

      subroutine b2agx0 (nrid, nCi, nFc, nVx)
      use b2mod_types
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nrid
*   ..output arguments (unspecified on entry)
      integer nCi, nFc, nVx

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2AGX0 reads the mesh dimensions from a formatted file.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine reads from unit nrid the mesh dimensions, nCi, nFc and
*     nVx.
*     Unit nrid must be connected to a formatted file. The input file
*     must be positioned properly: the first record must begin with the
*     text '*format', followed by 'unstructured'. The next record 
*     must begin with '*dimens ' or '*DIMENS ', and must contain the 
*     values of nCi, nFc and nVx in free format. Only those three
*     records are read.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nrid - integer, input.
*     nrid specifies a unit number for formatted input.
*
*     nCi, nFc, nVx - integer, output.
*     nCi, nFc and nVx specify the number of interior cells, the 
*     number of faces and the number of vertices respectively.   
*
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
      call subini ('b2agx0')
*   ..test nrid
      call xertst (1.le.nrid, 'faulty argument nrid')
*   ..read header line
      read (nrid,'(a8)',err=91) id
      call xertst (streql(id,'*dimens'),
     &  'faulty input id at *dimens, id='//id)
*   ..read nCi, nFc, nVx
      read (nrid,*,err=92) nCi, nFc, nVx
      call xertst (1.le.nCi.and.1.le.nFc.and.1.le.nVx,
     &  'faulty input nCi, nFc, nVx')

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.error conditions

   91 call xerrab ('error trying to read header *dimens')
   92 call xerrab ('error trying to read after '//id)

*-----------------------------------------------------------------------
*.end b2agx0

      end subroutine b2agx0
