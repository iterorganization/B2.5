#ifdef JET
      logical function my_isnan(a)
      implicit none
      real*8 :: a
      if(a.ne.a) then
       my_isnan=.true.
      else
       my_isnan=.false.
      endif
      return
      end function
#endif

      subroutine nancheck(nx, ny, data, label)
      use b2mod_types
#ifndef JET
      use ieee_arithmetic    ! IGNORE
#endif
      implicit none
#ifdef JET
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
#ifndef JET
        if(ieee_is_nan(data(ix,iy))) then
#else
        if(my_isnan(data(ix,iy))) then
#endif
           write(*,*) label, ' NaN found at ', ix,iy
        endif
       enddo
      enddo

      return
      end

