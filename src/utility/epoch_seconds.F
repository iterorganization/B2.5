      function epoch_seconds()
      use b2mod_types
      implicit none
      real (kind=R8) :: epoch_seconds

      character*8 date
      character*10 time
      character*5 zone
      integer timval(8)
      integer DaysToMonth(13),iday

      data DaysToMonth /0,31,59,90,120,151,181,212,243,273,304,334,365/

      call date_and_time(date, time, zone, timval)

      iday=365 * (timval(1) - 1970) + DaysToMonth(timval(2))
     & + timval(3) -1 + (timval(1) - 1969) / 4
      if ((timval(2).gt.2).and.(mod(timval(1),4).eq.0)) iday=iday+1
      epoch_seconds = timval(8)/1000.0_R8 + timval(7) +
     1 60 * (timval(6)-timval(4)) + 3600 * (timval(5) + 24 * iday)

      return
      end function epoch_seconds
