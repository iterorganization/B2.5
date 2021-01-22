      subroutine import(impaux)
c
c  version : 28.05.95 18:19
c
#include "impcom.inc"
      character*3 w,bfin*8
      data bfin /'finish'/
c
      call bsend(bfin,fin,8)
      promp=repeat(' ',16)
      write (w,700) maxlin
      fr(2:4)=w
      write (w,700) lenw
      fd(6:8)=w
      fe(6:8)=w
      fi(6:8)=w
  700 format (i3)
c
      erenc=.false.
      ninp=0
      kline=0
c
   10 call implin(nam)
      if(nonam) go to 30
   20 if(nam.eq.fin) go to 40
      done=.false.
      call impaux(nam)
      erenc=erenc.or.dumanl.or.(.not.done)
      if(done) go to 25
      call impmes(6)
      yline = .false.
      go to 10
   25 if(yline) go to 20
      go to 10
c
   30 call impmes(2)
      erenc=.true.
      go to 10
c
   40 call impmes(0)
c
      return
      end subroutine import
c
c
c=====================================================
      block data impblk
c
#include "impcom.inc"
c
      data  jmess ,  mimp ,  mess , keylen,  lfin /
     /        2   ,   5   ,   6   ,   6   ,   72  /,
     ,       lena ,  mfr  , immed , fixkey/
     /        6   ,   1   ,.true. ,.false./
      data aform/'a','d','e','i'/
      data number/'1','2','3','4','5','6','7','8',
     ,   '9','0',' ',',','.','+','-','(',')'/
      data term / ' ', ',', '*', '(', ')' /
      data        fd     ,     fe       ,    fi      /
     /      '(bn,d   .0)', '(bn,e   .0)', '(bn,i   )'/,
     ,             fr        /
     /        '(   a1)'/
      data  len / 8,8,4,4 /
c=====================================================
      end block data impblk
