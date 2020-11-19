*-----------------------------------------------------------------------
*.specification

      subroutine b2agx0 (nrid, nx, ny, nx1, ny1)
      use b2mod_types
      implicit none
*   ..input arguments (unchanged on exit)
      integer nrid
*   ..output arguments (unspecified on entry)
      integer nx, ny, nx1, ny1

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
*     This routine reads from unit nrid the mesh dimensions, nx and ny,
*     and the auxiliary mesh dimensions, nx1 and ny1.
*     Unit nrid must be connected to a formatted file. The input file
*     must be positioned properly: the first record must begin with the
*     text '*dimens ' or '*DIMENS ', and the next record must contain
*     the values of nx and ny in free format. Only those two records
*     are read.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nrid - integer, input.
*     nrid specifies a unit number for formatted input.
*
*     nx, ny - integer, output.
*     nx and ny specify the number of interior cells along the first
*     and the second coordinate, respectively. The total number of
*     cells is (nx+2)*(ny+2); they are indexed by (-1:nx,-1:ny).
*     It will hold that 0.le.nx and 0.le.ny.
*
*     nx1, ny1 - integer.
*     nx1 and ny1 specify the mesh dimensions for an auxiliary mesh
*     that will be used to evaluate accurately the metric coefficients.
*     This auxiliary mesh has no role outside the present program.
*     nx1 will be a multiple of nx and ny1 will be a multiple of ny.
*     The multiplicity test is put off till b2agmx, so we know if there
*     are any internal isolating cuts inside the domain.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      character id*8
*   ..procedures
      logical streql
      external subini, subend, xertst, xerrab, streql

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
*   ..read nx, ny, nx1, ny1
      read (nrid,*,err=92) nx, ny, nx1, ny1
      call xertst (1.le.nx.and.1.le.ny, 'faulty input nx, ny')
      call xertst (1.le.nx1.and.1.le.ny1, 'faulty input nx1, ny1')

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
