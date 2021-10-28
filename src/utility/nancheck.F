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

      subroutine nancheck(nx, ny, data, label)
      use b2mod_types
#ifndef NO_IEEE
      use ieee_arithmetic    ! IGNORE
#endif
      implicit none
#ifdef NO_IEEE
      logical, external :: my_isnan
#endif
*   ..input arguments (unchanged on exit)
      integer nx, ny
      real (kind=R8) ::
     *   data(-1:nx,-1:ny)
      character (len=*) label
*   ..local variables
      integer ix, iy

      do iy = -1, ny
        do ix = -1, nx
#ifndef NO_IEEE
          if(ieee_is_nan(data(ix,iy))) then
#else
          if(my_isnan(data(ix,iy))) then
#endif
            write(*,*) label, ' NaN found at ', ix,iy
          endif
        enddo
      enddo

      return
      end subroutine nancheck

