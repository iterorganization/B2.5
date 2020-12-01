      subroutine b2ruzd(nget, version, ns, zamin, zamax, zn, am, todo)
      use b2mod_types
      use b2mod_b2cmpa
     , ,only : is_neutral, partition_number, nspecies, amtol
      implicit none
*   ..input arguments (unchanged on exit)
      integer nget, ns
      character*10 version
      logical todo
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  zamin(0:ns-1), zamax(0:ns-1), zn(0:ns-1), am(0:ns-1)

*   ..local variables
      integer is
      real (kind=R8) ::
     *  za(0:ns-1), za2(0:ns-1)

*   ..procedures
      external subini, subend, xertst, cfrure

*   ..main computation
      call subini ('b2ruzd')
*   ..test input
!pb      call xertst (1.le.nget, 'faulty input nget')
      call xertst (1.le.ns, 'faulty input ns')

      if (version.lt.'01.001.024') then
        call cfrure (nget, ns, za, 'za')
        call cfrure (nget, ns, zn, 'zn')
        call cfrure (nget, ns, am, 'am')
        call cfrure (nget, ns, za2, 'za2')
        do is = 0, ns-1
          zamin(is) = za(is)
          zamax(is) = za(is)
        end do
      else
        call cfrure (nget, ns, zamin, 'zamin')
        call cfrure (nget, ns, zamax, 'zamax')
        call cfrure (nget, ns, zn, 'zn')
        call cfrure (nget, ns, am, 'am')
      end if
      if (todo) then
        do is = 0, ns-1
          is_neutral(is) = zamin(is).eq.0.0_R8
        end do
        nspecies = 1
        is = 0
        partition_number(is) = 0
        do while (is.lt.ns-1)
          if (lnext(is,is+1)) then
            partition_number(is+1) = partition_number(is)+1
          else
            partition_number(is+1) = 0
            nspecies = nspecies + 1
          endif
          is = is + 1
        end do
      end if

      call subend ()
      return

      contains

      logical function lnext(is,js)
      implicit none
      integer, intent(in) :: is, js
      lnext = zamax(is).lt.zamin(js).and.zn(is).eq.zn(js).and.
     &  abs(am(is)-am(js)).lt.amtol
      return
      end function lnext

      end subroutine b2ruzd
