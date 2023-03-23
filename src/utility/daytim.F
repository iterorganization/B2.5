      character*17 function daytim ()
      use b2mod_types
      implicit none
*     ------------------------------------------------------------------
*     DAYTIM returns the date and time in a character*17 format.
*     This version employs uses the date_and_time intrinsic function
*     and returns a date of the form 2015-07-09  15:05
*     ------------------------------------------------------------------
      character*8 date
      character*10 time
      character*5 zone
      integer values(8)
*     ------------------------------------------------------------------
      call date_and_time (date, time, zone, values)
      write(daytim,'(i4.4,a1,i2.2,a1,i2.2,a2,i2.2,a1,i2.2)') 
     &     values(1), '-', values(2), '-', values(3), '  ',
     &     values(5), ':', values(6)
      return
*     ------------------------------------------------------------------
      end function daytim
