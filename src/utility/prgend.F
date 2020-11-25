      subroutine prgend ()
      use b2mod_types
      use b2mod_ipmain
      implicit none
*     ------------------------------------------------------------------
*     PRGEND performs any tasks required at the end of program
*     execution.
*     ------------------------------------------------------------------
      integer lwmx, lwimx
      external sysend, subend, lwmax, lwimax
*     ------------------------------------------------------------------
*   ..write maximum workspace use
*     (activate or deactivate as desired)
      if (.true.) then
       call lwmax (lwmx)
       call lwimax (lwimx)
       write (*,'(1x,a,i9)') 'maximum use of stackw%w: ', lwmx
       write (*,'(1x,a,i9)') 'maximum use of stackw%wi:', lwimx
      endif
*   ..test use of internal parameters
      call ipfini ()
*   ..subprogram shut-down calls
      call subend ()
*   ..system shut-down calls
      call sysend ()
#ifdef PERFMON
      call perfoff
      call perfout ('whole')
#endif
      return
*     ------------------------------------------------------------------
      end subroutine prgend
      subroutine prgfin ()
      use b2mod_types
      use b2mod_ipmain
      implicit none
*     ------------------------------------------------------------------
*     PRGFIN performs any tasks required at the end of program
*     execution.
*     ------------------------------------------------------------------
      integer lwmx, lwimx
      external sysend, subend, lwmax, lwimax
*     ------------------------------------------------------------------
*   ..write maximum workspace use
*     (activate or deactivate as desired)
      if (.true.) then
       call lwmax (lwmx)
       call lwimax (lwimx)
       write (*,'(1x,a,i9)') 'maximum use of stackw%w: ', lwmx
       write (*,'(1x,a,i9)') 'maximum use of stackw%wi:', lwimx
      endif
*   ..test use of internal parameters
      if (.false.) call ipfini ()
*   ..subprogram shut-down calls
      call subend ()
*   ..system shut-down calls
      call sysend ()
#ifdef PERFMON
      call perfoff
      call perfout ('whole')
#endif
      return
*     ------------------------------------------------------------------
      end subroutine prgfin
