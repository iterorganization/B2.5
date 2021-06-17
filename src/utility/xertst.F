      subroutine xertst (lcond, messg)
      use b2mod_types
      implicit none
      logical lcond
      character messg*(*)
*     ------------------------------------------------------------------
*     XERTST tests a condition; if violated, calls xerrab to terminate
*     execution.
*     ------------------------------------------------------------------
      external xerrab
*     ------------------------------------------------------------------
      if (.not.lcond) then
       call xerrab (messg)
      endif
      return
*     ------------------------------------------------------------------
      end

      subroutine b25_mpi_stop
#ifdef USE_MPI
#if MPI_VERSION > 2
      use mpi     ! IGNORE
#endif
#endif
      implicit none
#ifdef USE_MPI
#if MPI_VERSION < 3
      include 'mpif.h'
#endif
      integer ierr_mpi
#ifdef B25_EIRENE
      logical ltrigger

      ltrigger=.false.
      call mpi_bcast(ltrigger,1,MPI_LOGICAL,0,
     .                   MPI_COMM_WORLD,ierr_mpi)
#endif
      call mpi_barrier(MPI_COMM_WORLD,ierr_mpi)
      call mpi_finalize(ierr_mpi)
#endif
      return
      end subroutine b25_mpi_stop
