*-----------------------------------------------------------------------
*.specification

      program b2ab
      use b2mod_types
      use b2mod_xerset
      use b2mod_create_afn_bc_tables
      implicit none

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2AB (main program) prepares databases which allow to better capture
*     the effect of boundary conditions of the kinetic neutral models
*     (specifically, the TRIM database used by EIRENE) in fluid neutral models.
*     This program calls the driver routine b2abdr for the actual
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
*     The program b2ab prepares a database for a certain combination of
*     incident particle and wall material.
*     This program unit opens and closes the input/output units, and
*     may perform some other system-dependent operations. It calls
*     routine b2abdr to control the actual computation.
*
*     The input units are:
*     ninp(0): formatted; provides parameters describing the database.

*     The output units are:
*     nout(0): un*formatted; provides the database
*
*     (See routine b2cdca for the meaning of 'un*formatted'.)
*
*
*  4. references
*
*     /1/Horsten, Niels, Giovanni Samaey, and M. Baelmans.
*     Development and assessment of 2D fluid neutral models that include atomic databases and a microscopic reflection model.
*     Nuclear Fusion 57.11 (2017): 116043.
*     /2/Horsten, Niels, PhD Thesis, KU Leuven, 2019, Chapter 4
*
*     Important note: the advanced fluid neutral boundary conditions were originally
*     only intended for the hydrogenic atoms. This is also assumed in the cited references.
*     However, the same principle could perhaps also yield a significant improvement for the
*     recycling boundary conditions of impurity neutrals. This is one of the reasons the source code
*     to create the tables needed for the advanced fluid neutral boundary conditions are added to the
*     B2.5 source code.
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
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks

*   ..local variables
      integer ninpab(1), noutab(1:2), io
      logical file_open
*   ..procedures
      external prgini, prgend, cfopen
      external b2abdr
*   ..initialise input/output units
      data ninpab/50/
      data noutab/60,61/

*-----------------------------------------------------------------------
*.documentation-internal
*
*     The following common blocks have their outermost declaration in
*     this routine; they need not be preserved between calls.
*
*     Description of some local variables:
*
*     ninpab - (0) integer array.
*     ninpab specifies the input unit numbers.
*
*     noutab - (0:1) integer array.
*     noutab specifies the output unit numbers.
*
*
*-----------------------------------------------------------------------
*.computation

*   ..program start-up calls
      call prgini ('b2ab')
*   ..open files
      call cfopen (ninpab(1),'b2ab.dat','old','formatted')
      call cfopen (noutab(1),'b2ab.prt','new','formatted')
      call cfopen (noutab(2),'b2_afn_bcs','new','un*formatted')
*   ..mark output unit for error messages
      call xerset (1)
*   ..call principal driver routine
      call b2abdr (ninpab, noutab)
*   ..close other files
      do io=1,size(ninpab)
         inquire(unit=ninpab(io),opened=file_open)
         if (file_open) close(ninpab(io))
      enddo
      do io=1,size(noutab)
         inquire(unit=noutab(io),opened=file_open)
         if (file_open) close(noutab(io))
      enddo

* ..end of computation
      call prgend ()
      stop 'b2ab'

*-----------------------------------------------------------------------
*.end b2ab

      end program b2ab
