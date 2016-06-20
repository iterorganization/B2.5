      subroutine cfrure (nget, n, refun, id)
      use b2mod_types
      implicit none
      integer nget, n
      real (kind=R8) ::
     *   refun(0:n-1)
      character id*(*)
csw 04dec2014 added fmt data (as in cfwure) to avoid problems with
c             pgf90 compiler , see below
      integer ifmt
      character fmt(0:5)*32
      save fmt
      data fmt /'(1p,5(1x,e25.17e3))', '(1p,6(1x,e21.13e3))',
     &  '(1p,7(1x,e17.9e3))', '(1p,8(1x,e15.7e3))',
     &  '(1p,9(1x,e13.5e3))', '(1p,10(1x,e12.4e3))'/
*     ------------------------------------------------------------------
*     CFRURE reads a real*8 array from an un*formatted file.
*     ------------------------------------------------------------------
      integer n0, i, ios, iloc, jloc, ndigits
      character chf*12, idcod*8, idtyp*8, id0*32, id1*32, zeile*80
      logical streql
      external xertst, xerrab, streql
*     ------------------------------------------------------------------
      id0 = id
      inquire (nget,form=chf,err=91)
      if (streql(chf,'formatted')) then
*    ..formatted read
       read (nget,'(2a8,i12,4x,a32)',err=92) idcod, idtyp, n0, id1
       call xertst (streql('*cf:',idcod),
     &   'cfrure--not a header line, expected *cf:, found '//idcod)
       call xertst (streql(id0,id1),
     &   'cfrure--wrong label: expected '//id0//', found '//id1)
       call xertst (streql('real',idtyp),
     &   'cfrure--wrong type: expected real, found '//idtyp)
       call xertst (n.eq.n0, 'cfrure--dimension mismatch at '//id0)
       if (0.lt.n) then
        read (nget,'(a80)',err=93,iostat=ios) zeile
        backspace(nget)
        iloc = index(zeile,'.')
        jloc = index(zeile,'E')
        if (jloc.eq.0) jloc = index(zeile,'e')
        if (jloc.eq.0) jloc = index(zeile,'D')
        if (jloc.eq.0) jloc = index(zeile,'d')
        if (iloc.eq.0 .or. jloc.eq.0) goto 93
        ndigits = jloc-iloc-1
        select case (ndigits) 
        case (17)
          ifmt = 0
        case (13)
          ifmt = 1
        case (9)
          ifmt = 2
        case (7)
          ifmt = 3
        case (5)
          ifmt = 4
        case (4)
          ifmt = 5
        case default
          call xerrab ('Unrecognized data format, id='//id0)
        end select
csw 04dec2014 changed * --> fmt(ifmt)
        read (nget,fmt(ifmt),err=93,iostat=ios) (refun(i),i=0,n-1)
       endif
      else if (streql(chf,'unformatted')) then
*    ..unformatted read
       read (nget,err=94) idcod, idtyp, n0, id1
       call xertst (streql('*cf:',idcod),
     &   'cfrure--not a header line, expected *cf:, found '//idcod)
       call xertst (streql(id0,id1),
     &   'cfrure--wrong label: expected '//id0//', found '//id1)
       call xertst (streql('real',idtyp),
     &   'cfrure--wrong type: expected real, found '//idtyp)
       call xertst (n.eq.n0, 'cfrure--dimension mismatch at '//id0)
       if (0.lt.n) then
        read (nget,err=95) (refun(i),i=0,n-1)
       endif
      else
       call xerrab ('cfrure--invalid file form-property: '//chf)
      endif
      return
*     ------------------------------------------------------------------
   91 call xerrab ('cfrure--error return from inquire statement')
   92 call xerrab ('cfrure--read error at head, formatted, id='//id0)
   93 call xerrab ('cfrure--read error at data, formatted, id='//id0)
   94 call xerrab ('cfrure--read error at head, unformatted, id='//id0)
   95 call xerrab ('cfrure--read error at data, unformatted, id='//id0)
*     ------------------------------------------------------------------
      end subroutine cfrure

      subroutine cfrure_opt (nget, n, refun, id)
      use b2mod_types
      implicit none
      integer nget, n
      real (kind=R8) ::
     *   refun(0:n-1)
      character id*(*)
csw 04dec2014 added fmt data (as in cfwure) to avoid problems with
c             pgf90 compiler , see below
      integer ifmt
      character fmt(0:5)*32
      save fmt
      data fmt /'(1p,5(1x,e25.17e3))', '(1p,6(1x,e21.13e3))',
     &  '(1p,7(1x,e17.9e3))', '(1p,8(1x,e15.7e3))',
     &  '(1p,9(1x,e13.5e3))', '(1p,10(1x,e12.4e3))'/
*     ------------------------------------------------------------------
*     CFRURE_OPT reads a real*8 array from an un*formatted file.
*     The array may not be present, in which case the read is skipped.
*     If the array size is not matching the expected size, skips the
*     reading and returns a soft error.
*     ------------------------------------------------------------------
      integer n0, i, ios, iloc, jloc, ndigits
      character chf*12, idcod*8, idtyp*8, id0*32, id1*32, zeile*80
      logical streql
      external xertst, xerrab, streql
*     ------------------------------------------------------------------
      id0 = id
      inquire (nget,form=chf,err=91)
      if (streql(chf,'formatted')) then
*    ..formatted read
       read (nget,'(2a8,i12,4x,a32)',err=92,end=89) 
     &   idcod, idtyp, n0, id1
       call xertst (streql('*cf:',idcod),
     &   'cfrure_opt--not a header line, expected *cf:, found '//idcod)
       if (.not.streql(id0,id1)) then
         write(6,'(a)') 
     &    'In cfrure_opt: expected '//trim(id0)//', found '//trim(id1)
         write(6,*) 'Reading skipped!'
         backspace(nget)
         return
       endif
       call xertst (streql('real',idtyp),
     &   'cfrure_opt--wrong type: expected real, found '//idtyp)
       if (n.ne.n0) then
         write(6,'(a,i6,a,i6)')
     &    'In cfrure_opt: for array '//trim(id1)//', expected size ',n,
     &    ' but found size ',n0
         write(6,*) 'Reading skipped!'
         backspace(nget)
         return
       endif
       if (0.lt.n) then
        read (nget,'(a80)',err=93,iostat=ios) zeile
        backspace(nget)
        iloc = index(zeile,'.')
        jloc = index(zeile,'E')
        if (jloc.eq.0) jloc = index(zeile,'e')
        if (jloc.eq.0) jloc = index(zeile,'D')
        if (jloc.eq.0) jloc = index(zeile,'d')
        if (iloc.eq.0 .or. jloc.eq.0) goto 93
        ndigits = jloc-iloc-1
        select case (ndigits) 
        case (17)
          ifmt = 0
        case (13)
          ifmt = 1
        case (9)
          ifmt = 2
        case (7)
          ifmt = 3
        case (5)
          ifmt = 4
        case (4)
          ifmt = 5
        case default
          call xerrab ('Unrecognized data format, id='//id0)
        end select
        read (nget,fmt(ifmt),err=93,iostat=ios) (refun(i),i=0,n-1)
       endif
      else if (streql(chf,'unformatted')) then
*    ..unformatted read
       read (nget,err=94,end=89) idcod, idtyp, n0, id1
       call xertst (streql('*cf:',idcod),
     &   'cfrure_opt--not a header line, expected *cf:, found '//idcod)
       if (.not.streql(id0,id1)) then
         write(6,'(a)') 
     &    'In cfrure_opt: expected '//trim(id0)//', found '//trim(id1)
         write(6,*) 'Reading skipped!'
         backspace(nget)
         return
       endif
       call xertst (streql('real',idtyp),
     &   'cfrure_opt--wrong type: expected real, found '//idtyp)
       if (n.ne.n0) then
         write(6,'(a,i6,a,i6)')
     &    'In cfrure_opt: for array '//trim(id1)//', expected size ',n,
     &    ' but found size ',n0
         write(6,*) 'Reading skipped!'
         backspace(nget)
         return
       endif
       if (0.lt.n) then
        read (nget,err=95) (refun(i),i=0,n-1)
       endif
      else
       call xerrab ('cfrure_opt--invalid file form-property: '//chf)
      endif
      return
*     ------------------------------------------------------------------
   89 write(6,'(a)') 
     &    'In cfrure_opt: expected '//trim(id0)//', found end of file'
      write(6,*) 'Reading skipped!'
      return 
*     ------------------------------------------------------------------
   91 call xerrab ('cfrure_opt--error return from inquire statement')
   92 call xerrab ('cfrure_opt--read error at head, formatted, id='//
     .             id0)
   93 call xerrab ('cfrure_opt--read error at data, formatted, id='//
     .             id0)
   94 call xerrab ('cfrure_opt--read error at head, unformatted, id='//
     .             id0)
   95 call xerrab ('cfrure_opt--read error at data, unformatted, id='//
     .             id0)
*     ------------------------------------------------------------------
      end subroutine cfrure_opt


