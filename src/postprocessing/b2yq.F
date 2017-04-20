*-----------------------------------------------------------------------
*.specification

      program b2yq
      use b2mod_types
      implicit none
#ifdef USE_MPI
      include 'mpif.h'
#endif

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2YQ (main program) is a postprocessor for b2. It provides a
*     quick display of the evolution data. This program calls the driver
*     routine b2yqdr for the actual computation.
*
*
*  2. specification
*
*     Main program.
*
*
*  3. description (see also routine b2cdca)
*
*     The complete program b2yq provides a quick display of the
*     evolution of a b2 calculation.
*     This program unit opens and closes the input/output units, and
*     may perform some other system-dependent operations. It calls
*     routine b2yqdr to control the actual computation.
*
*     The input units are:
*     ninp(0): formatted; provides output control.
*     ninp(1): un*formatted; provides the evolution data.
*
*     The output units are:
*     nout(0): formatted; provides summary printed output.
*     nout(1): graphical; displays the evolution data.
*
*     A scratch unit is:
*     nscr(0): unformatted.
*
*     (See routine b2cdca for the meaning of 'un*formatted'.)
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
      integer ninp(0:1), nout(0:1), nscr(0:0), io
#ifdef USE_MPI
      integer rank_mpi,size_mpi,ierr_mpi
#endif
      logical file_open
*   ..procedures
      external prgini, prgend, xerset, cfopen, gnopen, gopks, gclks,
     &  gclwk
      external b2yqdr
*   ..initialise input/output units
      data ninp/50,51/
      data nout/60,61/
      data nscr/70/

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
*     nout - (0:1) integer array.
*     nout specifies the output unit numbers.
*
*     nscr - (0:0) integer array.
*     nscr specifies the scratch unit numbers.
*
*-----------------------------------------------------------------------
*.computation

*   ..program start-up calls
#ifdef USE_MPI
      call mpi_init(ierr_mpi)
      call mpi_comm_rank(MPI_COMM_WORLD,rank_mpi,ierr_mpi)
      call mpi_comm_size(MPI_COMM_WORLD,size_mpi,ierr_mpi)
      if(rank_mpi .ne. 0) then
        call mpi_barrier(MPI_COMM_WORLD,ierr_mpi)
        call mpi_finalize(ierr_mpi)
        stop 'b2yq_mpi_child'
      endif
#endif
      call prgini ('b2yq')
*   ..open gks
      call gopks (0, 0)
*   ..open files
      call cfopen (ninp(0),'b2yq.dat','old','formatted')
      call cfopen (ninp(1),'b2ftrace','old','un*formatted')
      call cfopen (nout(0),'b2yq.prt','replace','formatted')   !srv 29.04.09
      call gnopen (nout(1),'b2yq.plt')
      open (nscr(0),status='scratch',form='unformatted')
*   ..mark output unit for error messages
      call xerset (0)
*   ..call principal driver routine
      call b2yqdr (ninp, nout, nscr)
*   ..close graphics files
      call gclwk (nout(1))
*   ..close gks
      call gclks ()
*   ..close other files
      do io=0,size(ninp)-1
         inquire(unit=ninp(io),opened=file_open)
         if (file_open) close(ninp(io))
      enddo
      do io=0,size(nout)-1
         inquire(unit=nout(io),opened=file_open)
         if (file_open) close(nout(io))
      enddo
      do io=0,size(nscr)-1
         inquire(unit=nscr(io),opened=file_open)
         if (file_open) close(nscr(io))
      enddo

* ..end of computation
      call prgend ()
#ifdef USE_MPI
      call mpi_barrier(MPI_COMM_WORLD,ierr_mpi)
      call mpi_finalize(ierr_mpi)
#endif
      stop 'b2yq'

*-----------------------------------------------------------------------
*.end b2yq

      end program b2yq
