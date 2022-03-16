c
      subroutine symbol(textline,iun,ipos,ikind,ival,rval,sval,cval)
      use b2mod_types
      implicit none
      character*500 textline
      character*1 c
      integer iun, ipos, ikind,ival,getint
      real (kind=R8) ::
     *   rval
      character*8 cval
      character*(*) sval
      integer icnt, noofdigits
      logical negative
      logical isadigit,isplusorminus,isnegative,is_comment
      external isadigit,isplusorminus,isnegative,is_comment
    1 call skipblanks(textline,iun,ipos)
      ival=0
      rval=0.0_R8
      sval=' '
      cval=' '
      c=textline(ipos:ipos)
      if(is_comment(textline).or.c.eq.'!') then
        call readline(textline,iun,ipos)
        goto 1
      endif
      if(c.eq.'''' .or. c.eq.'"') then
        ikind=3
        ipos=ipos+1
        icnt=1
    2   do while(ipos.lt.500 .and. textline(ipos:ipos).ne.c)
          if(icnt.lt.len(sval)) sval(icnt:icnt)=textline(ipos:ipos)
          icnt=icnt+1
          ipos=ipos+1
        enddo
        if(textline(ipos:ipos).eq.c) then
          ipos=ipos+1
          if(ipos.le.500 .and. textline(ipos:ipos).eq.c) then
            if(icnt.lt.len(sval)) sval(icnt:icnt)=textline(ipos:ipos)
            icnt=icnt+1
            ipos=ipos+1
            goto 2
          endif
          return
        endif
        if(textline(ipos:ipos).ne.c .and. ipos.gt.500) then
          call readline(textline,iun,ipos)
          write(*,*)
     1     'String not terminated on line --- continuing on next line'
          goto 2
        endif
      endif
      if(isadigit(c) .or. isplusorminus(c) .or. c.eq.'.') then
c we have a number
        negative=isnegative(textline,iun,ipos)
        if(c.eq.'.') goto 5
        ikind=1
        ival=getint(textline,iun,ipos,noofdigits)
        if(ipos.eq.1) then
          if(negative) ival=-ival
          return
        endif
        c=textline(ipos:ipos)
    5   if(c.eq.'.' .or. c.eq.'e' .or. c.eq.'E') then
          ikind=2
          rval=ival
          if(c.eq.'e' .or. c.eq.'E') goto 10
          ipos=ipos+1
          if(ipos .gt. 500) call readline(textline,iun,ipos)
          ival=getint(textline,iun,ipos,noofdigits)
          rval=rval+real(ival)/10.0_R8**noofdigits
          if(negative) rval=-rval
          if(ipos .eq. 1) return
          c=textline(ipos:ipos)
   10     if(c.eq.'e' .or. c.eq.'E') then
            if(negative .and. rval .gt. 0.0_R8) rval=-rval
            ipos=ipos+1
            if(ipos .gt. 500) call readline(textline,iun,ipos)
            negative=isnegative(textline,iun,ipos)
            ival=getint(textline,iun,ipos,noofdigits)
            if(negative) ival=-ival
c	    write(*,*) rval,ival
            rval=rval*10.0_R8**real(ival)
c	    write(*,*) rval
          endif
        else
          if(negative) ival=-ival
        endif
        return
      else
        ikind=4
        cval=c
        icnt=1
   20   ipos=ipos+1
        if(ipos .gt. 500) call readline(textline,iun,ipos)
        c=textline(ipos:ipos)
        if(c.eq.' ') return
        icnt=icnt+1
        if(icnt.le.8) cval(icnt:icnt)=c
        goto 20
      endif
      end
