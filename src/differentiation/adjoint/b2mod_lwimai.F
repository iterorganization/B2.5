      module b2mod_lwimai
      use b2mod_types
      use b2mod_stack
      implicit none
      private
      public :: lwiset, lwiget, lwimax

*     ------------------------------------------------------------------
*     LWIMAI is a group of routines (entry points in lwimai) that
*     maintain a base address for free space in the workspace area
*     b2mod_stack%wi. The subroutines are lwiset, lwiget, and lwimax.
*     LWISET frees workspace by resetting the base address.
*     LWIGET reserves workspace by incrementing the base address.
*     LWIMAX returns the maximum workspace use.
*     ------------------------------------------------------------------
      integer lwrki, lwrkimx
      save lwrki, lwrkimx
      data lwrki/0/, lwrkimx/0/
*     Invariant: 0.le.lwrki.le.lwrkimx.le.nstcwi
*     ------------------------------------------------------------------
      contains
*     ------------------------------------------------------------------
      subroutine lwiset (lwi0)
      implicit none
      integer, intent(in) :: lwi0
      external prvrti, b25_mpi_stop
*     (revert base address to lwi0)
      if (lwi0.lt.0.or.lwrki.lt.lwi0) then
        write (0,*) 'lwiset--faulty argument lwi0'
        call b25_mpi_stop
      end if
*     (set newly-free space to 'not-a-number')
      call prvrti (lwrki-lwi0, wi(lwi0), 1)
      lwrki = lwi0
      return
      end subroutine lwiset
*     ------------------------------------------------------------------
      subroutine lwiget (nwi, lwi1)
      implicit none
      integer, intent(in) :: nwi
      integer, intent(out) :: lwi1
      intrinsic max
      external b25_mpi_stop
      if (nwi.lt.0) then
        write (0,*) 'lwiget--faulty argument nwi'
        call b25_mpi_stop
      end if
      lwi1 = lwrki
      lwrki = lwrki+nwi
      lwrkimx = max(lwrki,lwrkimx)
      if (nstcwi.lt.lwrki) then
        write (0,*) 'lwiget--insufficient workspace stackw%wi(0:*)'
        call b25_mpi_stop
      end if
      return
      end subroutine lwiget
*     ------------------------------------------------------------------
      subroutine lwimax (lwi1)
      implicit none
      integer, intent(out) :: lwi1
      lwi1 = lwrkimx
      return
      end subroutine lwimax
*     ------------------------------------------------------------------
      end module b2mod_lwimai
