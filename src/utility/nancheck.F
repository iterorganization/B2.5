      subroutine nancheck(nx, ny, data, label)
      use b2mod_types
      use ieee_arithmetic    ! IGNORE
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny
      real (kind=R8) ::
     *   data(-1:nx,-1:ny)
      character (len=*) label
*   ..local variables
      integer ix, iy

      do iy = -1, ny
       do ix = -1, nx
        if(ieee_is_nan(data(ix,iy))) then
           write(*,*) label, ' NaN found at ', ix,iy
        endif
       enddo
      enddo

      return
      end
