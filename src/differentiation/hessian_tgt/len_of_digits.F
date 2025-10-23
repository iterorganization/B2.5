      function len_of_digits(i)
      implicit none
      integer i, j
      integer len_of_digits

      len_of_digits=1
      j = i
      if(j.lt.0) then
        j = -j
        len_of_digits = len_of_digits+1
      end if
      do while (j.ge.10)
        j = j/10
        len_of_digits = len_of_digits+1
      end do
      return
      end function len_of_digits
