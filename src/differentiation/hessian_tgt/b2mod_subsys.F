      module b2mod_subsys
      use b2mod_lwimai
      use b2mod_lwmain
      use b2mod_openmp
      implicit none
*     ------------------------------------------------------------------
*     SUBSYS is a group of routines that
*     perform certain functions related to workspace management and
*     error handling. The subroutines are subini, subend, and subxnd.
*     SUBINI performs any tasks required at the start of subprogram
*     execution.
*     SUBEND performs any tasks required at the end of subprogram
*     execution.
*     SUBXND is employed to build down the stack when the program
*     is terminated because of an error condition.
*     ------------------------------------------------------------------
      integer nchmx, nchs, lcr0, lci0
      parameter (nchmx=100)
      character chstck(0:nchmx-1)*8
      integer lcstck(0:1,0:nchmx-1)
      save nchs, chstck, lcstck
      data nchs/0/
*     ------------------------------------------------------------------
*!!   Note: nchmx may be changed. nchmx=100 is more than enough to
*     handle any nonrecursive call chain that one might meet.
*     ------------------------------------------------------------------
      contains
*     ------------------------------------------------------------------
      subroutine subini (srname)
      implicit none
      character srname*(*)
      external xertst
#ifdef PERFMON
      call perfon (srname)
#endif
#ifdef _OPENMP
      if (get_num_threads().gt.1.or.in_parallel()) return !!FR do not do the stack tracing if running in parallel
#endif
*   ..test nchs
      call xertst (nchs.ne.nchmx, 'subini--error, nchs.eq.nchmx')
*   ..get base addresses for free workspace
      call lwget (0, lcr0)
      call lwiget (0, lci0)
*   ..push data onto the stack
      chstck(nchs) = srname
      lcstck(0,nchs) = lcr0
      lcstck(1,nchs) = lci0
      nchs = nchs+1
cxpb      write(*,'(a,a8,i3)') 'subini: ',trim(srname),nchs
      return
      end subroutine subini
*     ------------------------------------------------------------------
      subroutine subend ()
      implicit none
      external xertst

#ifdef _OPENMP
      if (get_num_threads().gt.1.or.in_parallel()) return  !!FR do not do the stack tracing if running in parallel
#endif
*   ..test nchs
      call xertst (nchs.gt.0, 'subend--error, nchs.le.0')
*   ..pop data from the stack
*     (lcr0 and lci0 are the base addresses at invocation of the present
*     subprogram)
      chstck(nchs) = ' '
      nchs = nchs-1
      lcr0 = lcstck(0,nchs)
      lci0 = lcstck(1,nchs)
cxpb      write(*,'(a,a8,i3)') 'subend: ',trim(chstck(nchs)),nchs
*   ..reset base addresses for free workspace
      call lwset (lcr0)
      call lwiset (lci0)
#ifdef PERFMON
      call perfoff
#endif
      return
      end subroutine subend
*     ------------------------------------------------------------------
      subroutine subxnd (lbot, srnam1)
      implicit none
      logical, intent(out) :: lbot
      character*8, intent(out) :: srnam1
*   ..build down the stack
      if (nchs.eq.0) then
       lbot = .true.
       srnam1 = ' '
      else
       nchs = nchs-1
       lbot = .false.
       srnam1 = chstck(nchs)
       lcr0 = lcstck(0,nchs)
       lci0 = lcstck(1,nchs)
       call lwset (lcr0)
       call lwiset (lci0)
      endif
      return
      end subroutine subxnd
*     ------------------------------------------------------------------
      end module b2mod_subsys
