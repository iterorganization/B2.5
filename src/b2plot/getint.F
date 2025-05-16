c
      integer function getint(textline,iun,ipos,noofdigits)
      use b2mod_types
      integer iun, ipos
      character*500 textline
      character*1 c
      integer ival, noofdigits
      logical isadigit
      c=textline(ipos:ipos)
      ival=0
      getint=ival
      if(.not.isadigit(c)) return
      noofdigits=1
      ival=ichar(c)-ichar('0')
   10 continue
      ipos=ipos+1
      if(ipos .gt. 500) goto 20
      c=textline(ipos:ipos)
      if(isadigit(c)) then
        noofdigits=noofdigits+1
        ival=ival*10+ichar(c)-ichar('0')
        goto 10
      endif
   20 continue
      getint=ival
      if(ipos .gt. 500) call readline(textline,iun,ipos)
      if(noofdigits .gt. 9) then
        write(*,*) 'Too many digits'
        write(*,'(a500)') textline
        textline=' '
        textline(ipos:ipos)='^'
        write(*,'(a500)') textline
        textline=' '
        if(ipos .ne. 1) call readline(textline,iun,ipos)
      endif
      return
      end function getint
