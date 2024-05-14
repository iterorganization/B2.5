      subroutine b2wuzd(nput, version, ns, zamin, zamax, zn, am)
      use b2mod_types
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nput, ns
      character*10 version 
      real (kind=R8) :: 
     *  zamin(0:ns-1), zamax(0:ns-1), zn(0:ns-1), am(0:ns-1)

*   ..local variables
      integer is
      real (kind=R8) ::
     *  za(0:ns-1), za2(0:ns-1)

*   ..procedures
      external xertst, cfwure

*   ..main computation
      call subini ('b2wuzd')
*   ..test input
!pb       call xertst (1.le.nput, 'faulty input nput')
      call xertst (1.le.ns, 'faulty input ns')
      
      if (version.lt.'01.001.024') then
        do is = 0, ns-1
          za(is) = (zamin(is)+zamax(is))/2.0_R8
          za2(is) = za(is)**2
        end do
        call cfwure (nput, ns, za, 'za')
        call cfwure (nput, ns, zn, 'zn')
        call cfwure (nput, ns, am, 'am')
        call cfwure (nput, ns, za2, 'za2')
      else
        call cfwure (nput, ns, zamin, 'zamin')
        call cfwure (nput, ns, zamax, 'zamax')
        call cfwure (nput, ns, zn, 'zn')
        call cfwure (nput, ns, am, 'am')
      end if

      call subend ()
      return
      end subroutine b2wuzd
