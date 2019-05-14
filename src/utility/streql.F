      logical function streql (str0, str1)
      implicit none
      character str0*(*), str1*(*)
*     ------------------------------------------------------------------
*     STREQL checks for equality between two strings, ignoring case
*     shift and trailing blanks.
*!!   The text of this routine is case-sensitive.
*     ------------------------------------------------------------------
      integer len0, len1, lenl, i, k
      character ch0, ch1
      intrinsic len, ichar
      external b25_mpi_stop
      logical upcase, nlcase, chreql, isnull
      upcase(ch0,ch1) =
     &  ichar('a').le.ichar(ch0).and.ichar(ch0).le.ichar('z').and.
     &  ichar(ch1)-ichar(ch0).eq.ichar('A')-ichar('a')
      nlcase(ch0,ch1) =
     & (ichar(ch0).eq.0.or.ichar(ch0).eq.ichar(' ')).and.
     & (ichar(ch1).eq.0.or.ichar(ch1).eq.ichar(' '))
      chreql(ch0,ch1) = ch0.eq.ch1.or.upcase(ch0,ch1).or.upcase(ch1,ch0)
     &                            .or.nlcase(ch0,ch1)
*     ------------------------------------------------------------------
*   ..test installation
*     (protection against errors of case conversion.)
      if (ichar('a').eq.ichar('A')) then
       write (0,*) 'streql--installation error'
       call b25_mpi_stop
      endif
*   ..preliminaries
      len0 = len(str0)
      len1 = len(str1)
*   ..search for mismatch
      k = 0
    1 continue
      if (k.lt.min(len0,len1)) then
       if (chreql(str0(k+1:k+1),str1(k+1:k+1))) then
        k = k+1
        goto 1
       endif
      endif
*   ..set return value
      isnull=.true.
      lenl = max(len0,len1)-min(len0,len1)
      if (len0.eq.len1) then
       streql = k.eq.len0
      else if (len0.lt.len1) then
       do i = 1, lenl
         isnull = isnull.and.ichar(str1(len0+i:len0+i)).eq.0
       enddo
       streql = k.eq.len0.and.(str1(len0+1:len1).eq.' '.or.isnull)
      else if (len1.lt.len0) then
       do i = 1, lenl
         isnull = isnull.and.ichar(str0(len1+i:len1+i)).eq.0
       enddo
       streql = k.eq.len1.and.(str0(len1+1:len0).eq.' '.or.isnull)
      endif
      return
*     ------------------------------------------------------------------
      end
