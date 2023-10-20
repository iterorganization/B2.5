      logical function strmas (str0, str1)
      use b2mod_types
      implicit none
      character str0*(*), str1*(*)
*     ------------------------------------------------------------------
*     STRMAS checks for a certain matching between two strings.
*     strmas(str0,str1) returns the truth-value of
*       (str1 matches a substring of str0) .
*     Here, two strings match if they are of equal length and match
*     character for character; two characters match if they are equal,
*     or if they are alphabetic and differ only in case, or if either
*     one is the wildcard character '.'.
*!!   The text of this routine is case-sensitive.
*     ------------------------------------------------------------------
      integer len0, len1, i, k
      intrinsic len, ichar
      external xertst
*     ------------------------------------------------------------------
*   ..test installation
*     (protection against errors of case conversion.)
      call xertst (ichar('a').ne.ichar('A'),
     . 'strmas--installation error')
*   ..preliminaries
      len0 = len(str0)
      len1 = len(str1)
*   ..search for match
      i = 0
      k = 0
*     (loop invariant: str1(1:k) matches str0(i+1:i+k), there is no
*     earlier match of length len1.)
    1 continue
      if (i.le.len0-len1.and.k.lt.len1) then
       if (cmatch(str0(i+k+1:i+k+1),str1(k+1:k+1))) then
        k = k+1
       else
        i = i+1
        k = 0
       endif
       goto 1
      endif
*   ..set return value
      strmas = k.eq.len1
      return
*     ------------------------------------------------------------------

      contains

      logical function upcase(ch0,ch1)
      implicit none
      character ch0, ch1
      intrinsic ichar
      upcase =
     &  ichar('a').le.ichar(ch0).and.ichar(ch0).le.ichar('z').and.
     &  ichar(ch1)-ichar(ch0).eq.ichar('A')-ichar('a')
      return
      end function upcase

      logical function cmatch(ch0,ch1)
      implicit none
      character chwc, ch0, ch1
      parameter (chwc='.')

      cmatch = ch0.eq.ch1.or.ch0.eq.chwc.or.ch1.eq.chwc.or.
     &  upcase(ch0,ch1).or.upcase(ch1,ch0)
      return
      end function cmatch

      end function strmas
