      subroutine implin(nme)
c
c  version : 12.02.93 19:06
c
#include "impcom.inc"
      character nme(maxkey),bbl*8
      data bbl / '        ' /
c=====================================================
c
      if(promp(1:8).ne.bbl .or. promp(9:16).ne.bbl)
     -   write (mess,700) promp(1:8),promp(9:16)
  700 format(1x,2a8)
      read (mimp,fr,err=40,end=35) (ln(i),i=1,lfin)
      yline = .true.
      kline=kline+1
      nonam=.true.
      lbeg=1
      do i=1,nnum
        if(ln(1).eq.number(i)) return
      end do
c
c*** a new name
c
      call bfill(nam,maxkey,blnc)
      if(.not.fixkey) then
c-----------------------------------------------------
c*** free key format
c
        do i=1,maxkey
          if(ln(i).eq.blnc) exit
          nme(lbeg)=ln(i)
          lbeg=lbeg+1
        end do
      else    
c-----------------------------------------------------
c*** fixed key format
c
        call bsend(ln,nam,keylen)
        lbeg=keylen+1
      end if
      nonam=.false.
      return
c=====================================================
   35 backspace mimp
   40 nam=fin
      nonam=.false.
      return
c=====================================================
      end subroutine implin
