c
      subroutine species(is,string,decorate)
      use b2mod_types
      use b2mod_b2plot
      use b2mod_elements
      use b2mod_b2cmpa
      implicit none
      character*(*) string
      integer is,izn,izmax,izmin,imi,len,lnblnk
      logical decorate
      external lnblnk

      izn=nint(zn(is))
      izmin=nint(zamin(is))
      izmax=nint(zamax(is))
      imi=nint(am(is))
      if(izn .eq. 0) izn=1
      if(izn .eq. 1) then
        if(imi .eq. 1) string=' H'
        if(imi .eq. 2) string=' D'
        if(imi .eq. 3) string=' T'
      else
        string=' ' // elements(izn)
      endif
      len=lnblnk(string)
      if(decorate) string(len+1:len+2)='^{'
      len=lnblnk(string)
      if(izmin .gt. 0) string(len+1:len+1)='+'
      len=lnblnk(string)
      if(izmin .lt. 10) then
        write(string(len+1:len+1),'(i1)')
     1   izmin
      else
        write(string(len+1:len+2),'(i2)')
     1   izmin
      endif
      len=lnblnk(string)
      if(izmax .gt. izmin) then
        string(len+1:len+1)='_'
        if(izmax .gt. 0) string(len+2:len+2)='+'
        len=lnblnk(string)
        if(izmax .lt. 10) then
          write(string(len+1:len+1),'(i1)')
     1     izmax
        else
          write(string(len+1:len+2),'(i2)')
     1     izmax
        endif
      endif
      len=lnblnk(string)
      if(decorate) string(len+1:len+1)='}'
      return
      end subroutine species
