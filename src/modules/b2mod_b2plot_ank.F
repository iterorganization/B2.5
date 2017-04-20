      module b2mod_b2plot_ank
      use b2mod_types
      implicit none

      integer stacklen
      parameter (stacklen=1000)
      integer istack(stacklen), iistack, irstack, isstack
      integer stack_req(7)
      real (kind=R8) :: rstack(stacklen)
      character*256 sstack(stacklen)
      real (kind=R8), allocatable :: data(:)
      integer ndata

      contains

      integer function ipop(istack,iistack)
      implicit none
      integer istack(stacklen), iistack
      if(iistack .eq. 0) then
        write(*,*) 'istack underflow'
        ipop=0
      else
        ipop=istack(iistack)
        iistack=iistack-1
      endif
      return
      end function ipop

      subroutine ipush(istack,iistack,ival)
      implicit none
      integer istack(stacklen),ival,iistack
      if(iistack.eq.stacklen) then
        write(*,*) 'istack overflow'
      else
        iistack=iistack+1
        istack(iistack)=ival
      endif
      return
      end subroutine ipush
c
      function rpop(rstack,irstack)
      use b2mod_types
      implicit none
      real (kind=R8) :: rpop
      real (kind=R8) :: rstack(stacklen)
      integer irstack
      if(irstack .eq. 0) then
        write(*,*) 'rstack underflow'
        rpop=0.0_R8
      else
        rpop=rstack(irstack)
        irstack=irstack-1
      endif
      return
      end function rpop

      subroutine rpush(rstack,irstack,rval)
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   rstack(stacklen),rval
      integer irstack
      if(irstack.eq.stacklen) then
        write(*,*) 'rstack overflow'
      else
        irstack=irstack+1
        rstack(irstack)=rval
      endif
      return
      end subroutine rpush
c
      function spop(sstack,isstack)
      implicit none
      character*256 sstack(stacklen)
      character*256 spop
      integer isstack
      if(isstack .eq. 0) then
        write(*,*) 'sstack underflow'
        spop=' '
      else
        spop=sstack(isstack)
        isstack=isstack-1
      endif
      return
      end function spop

      subroutine spush(sstack,isstack,sval)
      use b2mod_types
      implicit none
      character*256 sstack(stacklen),sval
      integer isstack
      if(isstack.eq.stacklen) then
        write(*,*) 'sstack overflow'
      else
        isstack=isstack+1
        sstack(isstack)=sval
      endif
      return
      end subroutine spush

      subroutine pop_stack(stack_req)
      implicit none
      integer stack_req(7),is,itmp
      do is=1,stack_req(1)
        itmp=ipop(istack,iistack)
      enddo
      end subroutine pop_stack

      logical function bad_args(stack_req)
      use b2mod_b2plot
      implicit none
      integer stack_req(7)
      bad_args=.false.
      if(stack_req(1).gt.iistack) then
        write(*,*) 'Insufficient data on integer stack'
        bad_args= .true.
      endif
      if(stack_req(2).gt.irstack) then
        write(*,*) 'Insufficient data on real stack'
        bad_args= .true.
      endif
      if(stack_req(3).gt.isstack) then
        write(*,*) 'Insufficient data on string stack'
        bad_args= .true.
      endif
      if(stack_req(4).gt.iwrk) then
        write(*,*) 'Insufficient data on matrix stack'
        bad_args= .true.
      endif
      if(stack_req(5).gt.iewrk) then
        write(*,*) 'Insufficient data on Eirene stack'
        bad_args= .true.
      endif
      if(stack_req(6).gt.iwwrk) then
        write(*,*) 'Insufficient data on wall stack'
        bad_args= .true.
      endif
      if(stack_req(7).gt.itwrk) then
        write(*,*) 'Insufficient data on target stack'
        bad_args= .true.
      endif
      end function bad_args

      subroutine alloc_b2mod_b2plot_ank(nwall)
      implicit none
      integer nwall

      allocate(data(nwall))

      return
      end subroutine alloc_b2mod_b2plot_ank

      subroutine dealloc_b2mod_b2plot_ank
      implicit none

      deallocate(data)

      return
      end subroutine dealloc_b2mod_b2plot_ank

      end module b2mod_b2plot_ank
