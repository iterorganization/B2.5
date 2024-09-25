      subroutine prgini (prname)
      use b2mod_types
      use b2mod_stack
      use b2mod_subsys
      implicit none
      character prname*(*)
*     ------------------------------------------------------------------
*     PRGINI performs any tasks required at the start of program
*     execution.
*     ------------------------------------------------------------------
      external prvrt, prvrti
*     ------------------------------------------------------------------
*   ..system start-up calls
#ifdef PERFMON
      call perfinit
      call perfon ('whole')
#endif
      call sysini ()
*   ..subprogram start-up calls
      call subini (prname)
*   ..initialise w(0:nstckw-1) to 'not-a-number'
      call prvrt (nstckw, w, 1)
*   ..initialise wi(0:nstcwi-1) to 'not-a-number'
      call prvrti (nstcwi, wi, 1)
      return
*     ------------------------------------------------------------------
      end subroutine prgini
