#ifdef NO_IEEE
      logical function my_isnan(a)
      implicit none
      real*8 :: a
      if(a.ne.a) then
       my_isnan=.true.
      else
       my_isnan=.false.
      endif
      return
      end function my_isnan
#endif

      logical function my_isnanint(a)
      implicit none
      integer :: a
      if(a.ne.a) then
       my_isnanint=.true.
      else
       my_isnanint=.false.
      endif
      return
      end function my_isnanint

      subroutine nancheck(nCv, data, label)
      use b2mod_types
#ifndef NO_IEEE
      use ieee_arithmetic    ! IGNORE
#endif
      implicit none
#ifdef NO_IEEE
      logical, external :: my_isnan
#endif
*   ..input arguments (unchanged on exit)
      integer nCv
      real (kind=R8) :: data(nCv)
      character (len=*) label
*   ..local variables
      integer iCv

      do iCv = 1, nCv
#ifndef NO_IEEE
        if(ieee_is_nan(data(iCv))) then
#else
        if(my_isnan(data(iCv))) then
#endif
          write(*,*) label, ' NaN found at ', iCv
        endif
      enddo

      return
      end subroutine nancheck

      subroutine nancheck2d(nCv, istart, iend, data, label)
      use b2mod_types
#ifndef NO_IEEE
      use ieee_arithmetic    ! IGNORE
#endif
      implicit none
#ifdef NO_IEEE
      logical, external :: my_isnan
#endif
*   ..input arguments (unchanged on exit)
      integer nCv, istart, iend
      real (kind=R8) :: data(nCv,istart:iend)
      character (len=*) label
*   ..local variables
      integer i, iCv

      do i = istart, iend
        do iCv = 1, nCv
#ifndef NO_IEEE
          if(ieee_is_nan(data(iCv,i))) then
#else
          if(my_isnan(data(iCv,i))) then
#endif
            write(*,*) label, ' NaN found at ', iCv, i
          endif
        enddo
      enddo

      return
      end subroutine nancheck2d

      subroutine nancheckint(nCv, data, label)
      use b2mod_types
      implicit none
      logical, external :: my_isnanint
*   ..input arguments (unchanged on exit)
      integer nCv
      integer :: data(nCv)
      character (len=*) label
*   ..local variables
      integer iCv

      do iCv = 1, nCv
        if(my_isnanint(data(iCv))) then
          write(*,*) label, ' NaN found at ', iCv
        endif
      enddo

      return
      end subroutine nancheckint

      subroutine nancheck2dint(nCv, istart, iend, data, label)
      use b2mod_types
      implicit none
      logical, external :: my_isnanint
*   ..input arguments (unchanged on exit)
      integer nCv, istart, iend
      integer :: data(nCv,istart:iend)
      character (len=*) label
*   ..local variables
      integer i, iCv

      do i = istart, iend
        do iCv = 1, nCv
          if(my_isnanint(data(iCv,i))) then
            write(*,*) label, ' NaN found at ', iCv, i
          endif
        enddo
      enddo

      return
      end subroutine nancheck2dint
