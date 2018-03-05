*-----------------------------------------------------------------------
*.specification

      program b2ar
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
*     B2AR (main program) prepares a table of atomic rate coefficients
*     for B2. This program calls the driver routine b2ardr for the
*     actual computation.
*
*
*  2. specification
*
*     Main program.
*
*
*  3. description (see also routine b2cdca)
*
*     B2 is composed of several independent programs that communicate
*     through external files. b2ar is one of these programs; it
*     computes and writes to an external file a table of atomic rate
*     coefficients. The table may be computed either by use of the
*     ADPAK package developed by R. Hulse /3,4/, or by use of the
*     STRAHL package developed by K. Behringer /5,6/. Selected rate
*     coefficients for hydrogen will be re-evaluated using tabulated
*     data obtained from a code developed by J. C. Weisheit /7,8/.
*
*     This program unit opens and closes the input/output units, and
*     may perform some other system-dependent operations. It calls
*     routine b2ardr to control the actual computation.
*
*     The input units are:
*     ninp(0): formatted; provides species parameters.
*     ninp(1): formatted; provides input data for STRAHL calculation.
*     ninp(2): formatted; provides input data from Weisheit code.
*     ninp(3): un*formatted; provides default physics parameters.
*     ninp(4): un*formatted; provides the geometry.
*
*     The output units are:
*     nout(0): formatted; provides printed summary of rate coefs.
*     nout(1): un*formatted; provides the table of rate coefficients.
*
*     (See routine b2cdca for the meaning of 'un*formatted'.)
*
*
*  4. references
*
*     (Refs. /1/ and /2/ are out of date; !! to be replaced.)
*     /1/ B. J. Braams, Computational Studies in Tokamak Equilibrium
*     and Transport, thesis, State University at Utrecht, June 1986.
*     /2/ B. J. Braams, A Multi-Fluid Code for Simulation of the Edge
*     Plasma in Tokamaks, NET Report Nr 68, January 1987.
*     /3/ R. A. Hulse, Numerical Studies of Impurities in Fusion
*     Plasmas, Nuclear Technol./Fusion, vol.3 (1983) pp.259-272.
*     /4/ R. A. Hulse, !! (detailed description?).
*     /5/ K. Behringer, Description of the Impurity Transport Code
*     STRAHL, JET Report R(87)08, 1987.
*     /6/ K. Behringer, !! (journal article?).
*     /7/ Jon C. Weisheit, Recombination in Dense Plasmas, J. Phys. B
*     (Atom. Molec. Phys.) 8 (1975) 2556.
*     /8/ R. K. Janev, D. E. Post, W. D. Langer, K. Evans, D. B.
*     Heifetz and J. C. Weisheit, Survey of Atomic Processes in Edge
*     Plasmas, J. Nuclear Mater. 121 (1984) 10-16.
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
*
*  7. accuracy
*
*     Contact the author of the package used for a careful assessment.
*     The coefficients for ionisation and radiation should be accurate
*     (5-20 percent error?) at low density, but less so at high
*     density where multistage processes play a role. The coefficients
*     for recombination might be wrong by about 50 percent (?);
*     fortunately, volume recombination is often unimportant. The
*     coefficients for charge exchange may be off by 50 percent or
*     more (?).
*
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks

*   ..local variables
      integer ninp(0:4), nout(0:1)
      integer io
#ifdef USE_MPI
      integer rank_mpi,size_mpi,ierr_mpi
#endif
      logical file_open
*   ..procedures
      external prgini, prgend, xerset, cfopen
      external b2ardr
*   ..initialise input/output unit numbers
      data ninp/50,51,52,53,54/
      data nout/60,61/

*-----------------------------------------------------------------------
*.documentation-internal
*
*     The following common blocks have their outermost declaration in
*     this routine; they need not be preserved between calls.
*
*     Description of some local variables:
*
*     ninp - (0:4) integer array.
*     ninp specifies the input unit numbers.
*
*     nout - (0:1) integer array.
*     nout specifies the output unit numbers.
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
        stop 'b2ar_mpi_child'
      endif
#endif
      call prgini ('b2ar')
*   ..open files
      call cfopen (ninp(0),'b2ar.dat','old','formatted')
      call cfopen (ninp(1),'stra.dat','old','formatted')
      call cfopen (ninp(2),'weis.dat','old','formatted')
      call cfopen (ninp(3),'b2fpardf','old','un*formatted')
      call cfopen (ninp(4),'b2fgmtry','old','un*formatted')
      call cfopen (nout(0),'b2ar.prt','new','formatted')
      call cfopen (nout(1),'b2frates','new','un*formatted')
*   ..mark output unit for error messages
      call xerset (0)
*   ..compute, display and save table
      call b2ardr (ninp, nout)
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
      call prgend ()
#ifdef USE_MPI
      call mpi_barrier(MPI_COMM_WORLD,ierr_mpi)
      call mpi_finalize(ierr_mpi)
#endif
      stop 'b2ar'

*-----------------------------------------------------------------------
*.end b2ar

      end program b2ar
