      character*17 function daytim ()
      use b2mod_types
      implicit none
*     ------------------------------------------------------------------
*     DAYTIM returns the date and time in a character*17 format.
*!!   This routine may contain system dependencies.
*     This version employs NAG library routines and returns a result
*     of the form '199x Mon dd hh:mm'.
*     ------------------------------------------------------------------
      integer itime(0:6)
      character chtd*30, x05abf*30
      external x05aaf, x05abf
*     ------------------------------------------------------------------
      call x05aaf (itime)
      chtd = x05abf(itime)
      if (chtd(5:5).eq.' ') then
       chtd(5:5) = '0'
      endif
*     (chtd is of the form 'Day nnth Mon 199x hh:mm:ss.000')
      daytim = chtd(14:18)//chtd(10:13)//chtd(5:6)//chtd(18:23)
      return
*     ------------------------------------------------------------------
      end
