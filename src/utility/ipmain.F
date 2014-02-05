      subroutine ipmain (cnam, cval)
      use b2mod_types
      implicit none
      character cnam*(*), cval*(*)
*     ------------------------------------------------------------------
*     IPMAIN is a group of routines (entry points in ipmain) that
*     operate on a collection of names and values representing "internal
*     parameters" in a code. The entry points are ipsetc, ipgetc, and
*     ipfini.
*     IPSETC saves the value of an internal parameter.
*     IPGETC returns the value of an internal parameter.
*     IPFINI writes an error message for any parameter that has been
*     set but not read.
*     ------------------------------------------------------------------
*     In conformance with the Fortran spirit, case is not significant
*     in the names of the internal parameters. The argument cnam is
*     converted to lower case before being stored in the internal table
*     or being compared against other entries in the internal table.
*     In conflict with the Fortran spirit, leading and embedded blanks
*     are significant in the names of internal parameters, although
*     trailing blanks are not.
*     Any variable may be retrieved at most ipgtmx times, and may
*     not be set after it has been retrieved; these restrictions are
*     imposed to guard against inappropriate use. Because the routine
*     was not designed for efficiency, ipgtmx is set to a small integer
*     value. In any subroutine in which a certain parameter is
*     required, that parameter should be retrieved just once (on the
*     first call).
*     ------------------------------------------------------------------
*   ..local variables
      integer iptinc, iptbmx, ipgtmx, ipcnt, ip, i, k
      real (kind=R8) ::
     *   xdmx
      parameter (iptinc=233, iptbmx=377, ipgtmx=20, xdmx=0.8_R8)
      character ipnam(0:iptbmx-1)*80, ipval(0:iptbmx-1)*80, cnam0*80
      integer ipstat(0:iptbmx-1)
      save ipcnt, ipnam, ipval, ipstat
*   ..procedures
      intrinsic mod, ichar
      external xertst, chcase
*   ..initialisation
      data ipcnt/0/, ipstat/iptbmx*-1/
*     ------------------------------------------------------------------
*!!   Note: The parameter iptbmx defines the size of the internal
*     parameter table, and may be changed as needed.
*     The routine employs a hash table. The size of the table is given
*     by iptbmx. The parameter iptinc defines a simple hash function;
*     it is advisable to let iptinc and iptbmx be incommensurate, and
*     best if iptinc and iptbmx are successive Fibonacci numbers. It
*     is treated as an error if the table fill ratio exceeds the value
*     of the parameter xdmx (0.lt.xdmx.lt.1; suggested value 0.8).
*     The table is contained in the arrays ipnam, ipval, and ipstat.
*     (Forall ip : 0.le.ip.lt.iptbmx : -1.le.ipstat(ip).le.ipgtmx).
*     ipnam(ip) is the name of an internal parameter, ipval(ip) is its
*     value, and ipstat(ip) is the number of times that it has been
*     retrieved; ipstat(ip).eq.-1 if the parameter has not been set.
*     ipcnt is the number of parameters that have been set.
*     ------------------------------------------------------------------
      return
*     ------------------------------------------------------------------
      entry ipsetc (cnam, cval)
*   ..standardise the name
      cnam0 = cnam
      call chcase (0, cnam0)
*   ..evaluate the hash function
      k = 0
      do i = 1, len(cnam0)
       k = k+ichar(cnam0(i:i))
      enddo
      ip = mod(k*iptinc,iptbmx)
*   ..search in the table
    1 if (0.le.ipstat(ip)) then
       if (cnam0.ne.ipnam(ip)) then
        ip = mod(ip+1,iptbmx)
        goto 1
       endif
      endif
*   ..test presence
      if (0.le.ipstat(ip)) then
*    ..replace an entry
       call xertst (ipstat(ip).eq.0,
     &   'ipsetc--error, parameter has already been read: '//cnam0)
       ipval(ip) = cval
      else
*    ..extend the table
       call xertst (ipcnt.lt.xdmx*iptbmx,
     &   'ipsetc--error, internal parameter table overflow')
       ipnam(ip) = cnam0
       ipval(ip) = cval
       ipstat(ip) = 0
       ipcnt = ipcnt+1
      endif
      return
*     ------------------------------------------------------------------
      entry ipgetc (cnam, cval)
*   ..standardise the name
      cnam0 = cnam
      call chcase (0, cnam0)
*   ..evaluate the hash function
      k = 0
      do i = 1, len(cnam0)
       k = k+ichar(cnam0(i:i))
      enddo
      ip = mod(k*iptinc,iptbmx)
*   ..search in the table
    2 if (0.le.ipstat(ip)) then
       if (cnam0.ne.ipnam(ip)) then
        ip = mod(ip+1,iptbmx)
        goto 2
       endif
      endif
*   ..test presence
      if (0.le.ipstat(ip)) then
*    ..return the entry
       call xertst (ipstat(ip).lt.ipgtmx,
     &   'ipgetc--error, entry is read too often: '//cnam0)
       cval = ipval(ip)
       ipstat(ip) = ipstat(ip)+1
      else
*      (no action, return cval unchanged)
      endif
      return
*     ------------------------------------------------------------------
      entry ipfini ()
      do ip = 0, iptbmx-1
       if (ipstat(ip).eq.0) then
        write (*,'(1x,a,1x,a)')
     &    '*** IPMAIN--WARNING: parameter has been set but not read:',
     &    ipnam(ip)
       endif
      enddo
      return
*     ------------------------------------------------------------------
      end subroutine ipmain
