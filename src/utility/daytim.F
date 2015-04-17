      character*17 function daytim ()
      use b2mod_types
      implicit none
*     ------------------------------------------------------------------
*     DAYTIM returns the date and time in a character*17 format.
*     This version employs a Fortran intrinsic procedure and returns a
*     result of the form 'yyyy mm dd hh:mm.'.
*     ------------------------------------------------------------------
      integer values(8)
      character date*8, time*10, zone*5
*     ------------------------------------------------------------------
      call date_and_time(date,time,zone,values)
*     date has the form ccyymmdd
*     time has the form hhmmss.sss
      daytim = date(1:4)//' '//date(5:6)//' '//date(7:8)//' '
     $  //time(1:2)//':'//time(3:4)//'.'
*     daytim contains yyyy mm dd hh:mm
      return
*     ------------------------------------------------------------------
      end
