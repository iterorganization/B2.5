      module b2mod_lwmain
      use b2mod_types
      use b2mod_stack
      implicit none
      private
      public :: lwset, lwget, lwmax
*     ------------------------------------------------------------------
*     LWMAIN is a group of routines (entry points in lwmain) that
*     maintain a base address for free space in the workspace area
*     b2mod_stack%w. The subroutines are lwset, lwget, and lwmax.
*     LWSET frees workspace by resetting the base address.
*     LWGET reserves workspace by incrementing the base address.
*     LWMAX returns the maximum workspace use.
*     ------------------------------------------------------------------
      integer lwrk, lwrkmx
      save lwrk, lwrkmx
      data lwrk/0/, lwrkmx/0/
*     Invariant: 0.le.lwrk.le.lwrkmx.le.nstckw
*     ------------------------------------------------------------------
      contains
*     ------------------------------------------------------------------
      subroutine lwset (lw0)
      implicit none
      integer, intent(in) :: lw0
      external prvrt, b25_mpi_stop
*     (revert base address to lw0)
      if (lw0.lt.0.or.lwrk.lt.lw0) then
        write (0,*) 'lwset--faulty argument lw0'
        call b25_mpi_stop
      end if
*     (set newly-free space to 'not-a-number')
      call prvrt (lwrk-lw0, w(lw0), 1)
      lwrk = lw0
      return
      end subroutine lwset
*     ------------------------------------------------------------------
      subroutine lwget (nw, lw1)
      implicit none
      integer, intent(in) :: nw
      integer, intent(out) :: lw1
      intrinsic max
      external b25_mpi_stop
      if (nw.lt.0) then
        write (0,*) 'lwget--faulty argument nw'
        call b25_mpi_stop
      end if
      lw1 = lwrk
      lwrk = lwrk+nw
      lwrkmx = max(lwrk,lwrkmx)
      if (nstckw.lt.lwrk) then
        write(0,*) 'lwget--insufficient workspace b2mod_stack%w(0:*)'
        call b25_mpi_stop
      end if
      return
      end subroutine lwget
*     ------------------------------------------------------------------
      subroutine lwmax (lw1)
      implicit none
      integer, intent(out) :: lw1
      lw1 = lwrkmx
      return
      end subroutine lwmax
*     ------------------------------------------------------------------
      end module b2mod_lwmain
