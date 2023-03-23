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
      real (kind=R8) ::
     *   data(nCv)
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

