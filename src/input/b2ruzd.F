      subroutine b2ruzd(nget, version, ns, zamin, zamax, zn, am, todo)
      use b2mod_types
      use b2mod_b2cmpa
     , ,only : is_neutral, partition_number, nspecies
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nget, ns
      character*10 version
      logical todo
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  zamin(0:ns-1), zamax(0:ns-1), zn(0:ns-1), am(0:ns-1)

*   ..local variables
      integer is, js
      real (kind=R8) ::
     *  amtol, za(0:ns-1), za2(0:ns-1)
      parameter (amtol=1.0e-2_R8)

*   ..procedures
      intrinsic nint
      logical lnext
      external xertst, cfrure
      lnext(is,js) = zamax(is).lt.zamin(js).and.zn(is).eq.zn(js).and.
     &  abs(am(is)-am(js)).lt.amtol

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
      end subroutine b2ruzd
