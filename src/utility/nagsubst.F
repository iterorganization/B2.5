#ifdef NO_NAG
      subroutine e01bef (n, x, f, d, ifail)
      use b2mod_types
      implicit none
      integer n, ifail
      real (kind=R8) ::
     *   x(n), f(n), d(n)
      external subini, subend, xertst
*     ------------------------------------------------------------------
*     A crude substitute for the NAG routine of the same name, for the
*     case of an eight-byte real*8 environment.
*     ------------------------------------------------------------------
      integer i
      call subini ('e01bef')
      call xertst (n.ge.2, 'faulty argument n')
      do i = 1, n-1
        call xertst (x(i).lt.x(i+1), 'faulty argument x')
      enddo
      do i = 1, n-1
        d(i)=(f(i+1)-f(i))/(x(i+1)-x(i))
      enddo
      d(n)=0.0_R8
      ifail=0
      call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine e01bef
      subroutine e01bff (n, x, f, d, m, px, pf, ifail)
      use b2mod_types
      implicit none
      integer n, m, ifail
      real (kind=R8) ::
     *   x(n), f(n), d(n), px(m), pf(m)
      external subini, subend, xertst
*     ------------------------------------------------------------------
*     A crude substitute for the NAG routine of the same name, for the
*     case of an eight-byte real*8 environment.
*     ------------------------------------------------------------------
      integer i, j
      call subini ('e01bff')
      call xertst (n.ge.2, 'faulty argument n')
      do i = 1, n-1
        call xertst (x(i).lt.x(i+1), 'faulty argument x')
      enddo
      call xertst (m.ge.1, 'faulty argument m')
      ifail=0
      do j = 1, m
        if (px(j).lt.x(1)) then
          ifail = 4
          pf(j)=f(1)+d(1)*(px(j)-x(1))
        elseif (px(j).ge.x(n)) then
          ifail = 4
          pf(j) = f(n)+d(n-1)*(px(j)-x(n))
        else
          do i = 1, n-1
            if (px(j).ge.x(i) .and. px(j).lt.x(i+1)) then
              pf(j) = f(i)+d(i)*(px(j)-x(i))
            endif
          enddo
        endif
      enddo
      if (ifail.eq.4) then
        write(*,*) ' Extrapolation performed in e01bff !'
        write(*,*) ' Results unreliable !!!'
      endif
      call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine e01bff
      subroutine g02caf (n, x, y, result, ifail)
      use b2mod_types
      implicit none
      integer n, ifail
      real (kind=R8) ::
     *   x(n), y(n), result(20)
      integer i
*     ------------------------------------------------------------------
*     A crude substitute for the NAG routine of the same name, for the
*     case of an eight-byte real*8 environment.
*     ------------------------------------------------------------------
      ifail = 0
      result = 0.0_R8
      if (n.le.2) then
        ifail = 1
        return
      end if
      do i = 1, n
        result(1) = result(1) + x(1)
        result(2) = result(2) + y(1)
      end do
      result(1) = result(1)/real(n)
      result(2) = result(2)/real(n)
      do i = 1, n
        result(3) = result(3) + (x(i)-result(1))**2
        result(4) = result(4) + (y(i)-result(2))**2
      end do
      result(19) = result(4)
      do i = 1, n
        result(5) = result(5) + (x(i)-result(1))*(y(i)-result(2))
        result(6) = result(6) + (x(i)-result(1))*(y(i)-result(2))
      end do
      result(5) = result(5) / sqrt ( result(3)*result(4) )
      result(6) = result(6) / result(3)
      result(7) = result(2) - result(6)*result(1)
      do i = 1, n
        result(16) = result(16) + (y(i)-result(7)-result(6)*x(i))**2
      end do
      result(17) = n - 2
      result(18) = result(16) / result(17)
      result(20) = n - 1
      result(8) = sqrt ( result(18) / result(3) )
      result(9) = sqrt ( result(18) *
     .                 ( 1.0_R8/real(n) + result(1)**2/result(3) ) )
      result(10) = result(6) / result(8)
      result(11) = result(7) / result(9)
      result(12) = result(19) - result(16)
      result(13) = 1
      result(14) = result(12) / result(13)
      result(15) = result(14) / result(18)
      result(3) = sqrt ( result(3) / real(n-1) )
      result(4) = sqrt ( result(4) / real(n-1) )
      return
*     ------------------------------------------------------------------
      end subroutine g02caf
      subroutine g13abf (x, nx, nk, xm, xv, r, stat, ifail)
      use b2mod_types
      implicit none
      integer nx, nk, ifail
      real (kind=R8) ::
     *   x(nx), xm, xv, r(nk), stat
      integer i, k
*     ------------------------------------------------------------------
*     A crude substitute for the NAG routine of the same name, for the
*     case of an eight-byte real*8 environment.
*     ------------------------------------------------------------------
      ifail = 0
      xm = 0.0_R8
      xv = 0.0_R8
      stat = 0.0_R8
      r(1:nk) = 0.0_R8
      if (nx.le.nk .or. nx.le.1 .or. nk.le.0) then
        ifail = 1
        return
      end if
      do i = 1, nx
        xm = xm + x(i)
      end do
      xm = xm / real(nx)
      do i = 1, nx
        xv = xv + (x(i)-xm)**2
      end do
      do k = 1, nk
        do i = 1, nx - nk
          r(k) = r(k) + (x(i)-xm)*(x(i+k)-xm)
        end do
        r(k) = r(k) / xv
      end do
      xv = xv / real(nx-1)
      do k = 1, nk
        stat = stat + r(k)**2
      end do
      stat = nx*stat
      return
*     ------------------------------------------------------------------
      end subroutine g13abf
      function x02ajf ()
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   x02ajf
*     ------------------------------------------------------------------
*     A FORTRAN intrinsic for the NAG routine of the same name, for the
*     case of an eight-byte real*8 environment.
*     ------------------------------------------------------------------
      x02ajf = epsilon(1.0_R8)
      return
*     ------------------------------------------------------------------
      end function x02ajf
      function x02akf ()
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   x02akf
*     ------------------------------------------------------------------
*     A FORTRAN intrinsic for the NAG routine of the same name, for the
*     case of an eight-byte real*8 environment.
*     ------------------------------------------------------------------
      x02akf = tiny(1.0_R8)
      return
*     ------------------------------------------------------------------
      end function x02akf
      function x02alf ()
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   x02alf
*     ------------------------------------------------------------------
*     A FORTRAN intrinsic for the NAG routine of the same name, for the
*     case of an eight-byte real*8 environment.
*     ------------------------------------------------------------------
      x02alf = huge(1.0_R8)
      return
*     ------------------------------------------------------------------
      end function x02alf
      function x02amf ()
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   x02amf
*     --------------------------------------------------------------------
*     A FORTRAN intrinsic substitute for the NAG routine of the same name, 
*     for the case of an eight-byte real*8 environment.
*     --------------------------------------------------------------------
      x02amf = tiny(1.0_R8)
      return
*     ------------------------------------------------------------------
      end function x02amf
      subroutine x05aaf (itime)
      use b2mod_types
      implicit none
      integer itime(0:6)
*     --------------------------------------------------------------------
*     A FORTRAN intrinsic substitute for the NAG routine of the same name.
*     x05aaf should return the local date and time as a sequence of
*     integers, e.g. (1994, 01, 29, 16, 17, 18, 666) to represent
*     1994 Jan 29, 16:17:18.666 (24-hour clock to the millisecond).
*     --------------------------------------------------------------------
      character*8 date
      character*10 time
      character*5 zone
      integer values(8)
      call date_and_time (date, time, zone, values)
      itime(0) = values(1)
      itime(1) = values(2)
      itime(2) = values(3)
      itime(3) = values(5)
      itime(4) = values(6)
      itime(5) = values(7)
      itime(6) = values(8)
      return
*     ------------------------------------------------------------------
      end subroutine x05aaf
      character*30 function x05abf (itime)
      use b2mod_types
      implicit none
      integer itime(0:6)
*     --------------------------------------------------------------------
*     A FORTRAN intrinsic substitute for the NAG routine of the same name.
*     x05abf should return the date and time in the format
*     'Day ddth Mmm yyyy hh:mm:ss.000'
*     --------------------------------------------------------------------
      integer values(8)
      character*8 date
      character*10 time
      character*5 zone
      character*2 th
      character*3 month(12)
      data month /'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 
     .            'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/
      save month
      call date_and_time (date, time, zone, values)
      th = 'th'
      if (values(3).eq.1 .or. values(3).eq.21 .or. values(3).eq.31) 
     > th = 'st'
      if (values(3).eq.2 .or. values(3).eq.22) th = 'nd' 
      x05abf = 'Day '//date(7:8)//th//' '//month(values(2))//' '//
     . date(1:4)//' '//time(1:2)//':'//time(3:4)//':'//time(5:10)
      return
*     ------------------------------------------------------------------
      end function x05abf
#else
      subroutine REAL_NAG_ROUTINES_USED
      use b2mod_types
      end subroutine REAL_NAG_ROUTINES_USED
#endif
