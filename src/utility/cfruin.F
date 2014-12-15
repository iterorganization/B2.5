      subroutine cfruin (nget, n, infun, id)
      use b2mod_types
      implicit none
      integer nget, n, infun(0:n-1)
      character id*(*)
csw 04dec2014 added fmt data (as in cfwure) to acvoid problems with
c             pgf90 compiler , see below
      integer ifmt
      parameter (ifmt=2)
      character fmt(0:6)*32
      save fmt
      data fmt /'(8(1x,i15))', '(10(1x,i12))', '(12(1x,i10))',
     &  '(14(1x,i8))', '(16(1x,i7))', '(18(1x,i6))', '(22(1x,i5))'/
*     ------------------------------------------------------------------
*     CFRUIN reads an integer array from an un*formatted file.
*     ------------------------------------------------------------------
      integer n0, i
      character chf*12, idcod*8, idtyp*8, id0*32, id1*32
      logical streql
      external xertst, xerrab, streql
*     ------------------------------------------------------------------
      id0 = id
      inquire (nget,form=chf,err=91)
      if (streql(chf,'formatted')) then
*    ..formatted read
       read (nget,'(2a8,i12,4x,a32)',err=92) idcod, idtyp, n0, id1
       call xertst (streql('*cf:',idcod),
     &   'cfruin--not a header line, expected *cf:, found '//idcod)
       call xertst (streql(id0,id1),
     &   'cfruin--wrong label: expected '//id0//', found '//id1)
       call xertst (streql('int',idtyp),
     &   'cfruin--wrong type: expected int, found '//idtyp)
       call xertst (n.eq.n0, 'cfruin--dimension mismatch at '//id0)
       if (0.lt.n) then
csw 04dec2014 changed * --> fmt(ifmt)
        read (nget,fmt(ifmt),err=93) (infun(i),i=0,n-1)
       endif
      else if (streql(chf,'unformatted')) then
*    ..unformatted read
       read (nget,err=94) idcod, idtyp, n0, id1
       call xertst (streql('*cf:',idcod),
     &   'cfruin--not a header line, expected *cf:, found '//idcod)
       call xertst (streql(id0,id1),
     &   'cfruin--wrong label: expected '//id0//', found '//id1)
       call xertst (streql('int',idtyp),
     &   'cfruin--wrong type: expected int, found '//idtyp)
       call xertst (n.eq.n0, 'cfruin--dimension mismatch at '//id0)
       if (0.lt.n) then
        read (nget,err=95) (infun(i),i=0,n-1)
       endif
      else
       call xerrab ('cfruin--invalid file form-property: '//chf)
      endif
      return
*     ------------------------------------------------------------------
   91 call xerrab ('cfruin--error return from inquire statement')
   92 call xerrab ('cfruin--read error at head, formatted, id='//id0)
   93 call xerrab ('cfruin--read error at data, formatted, id='//id0)
   94 call xerrab ('cfruin--read error at head, unformatted, id='//id0)
   95 call xerrab ('cfruin--read error at data, unformatted, id='//id0)
*     ------------------------------------------------------------------
      end subroutine cfruin
