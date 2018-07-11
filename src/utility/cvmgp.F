      function cvmgp (x1,x2,x3)
      use b2mod_types
      implicit none
      real (kind=R8) :: cvmgp
      real (kind=R8) :: x1,x2,x3
*     ------------------------------------------------------------------
*     CVMGP( CRAY - routinen - ersatz).(RFS 27.5.1992)
*     x1 returned if x3 >= 0
*     x2 returned if x3 < 0
*     ------------------------------------------------------------------
      real (kind=R8) :: t
*     ------------------------------------------------------------------
      if (x3.ge.0.0_R8) then
       t = x1
      else
       t = x2
      endif
      cvmgp = t
      return
*     ------------------------------------------------------------------
      end function cvmgp
