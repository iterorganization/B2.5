c
      subroutine species(is,string,decorate)
      use b2mod_types
      use b2mod_b2plot
      use b2mod_elements
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
      character*(*) string
      integer is,izn,izmax,izmin,imi,len,lnblnk
      logical decorate
      external lnblnk

      call subini('species')
      izn=nint(zn(is))
      izmin=nint(zamin(is))
      izmax=nint(zamax(is))
      imi=nint(am(is))
      call xertst( izn.ge.0, 'Inconsistent nuclear charge input!')
      call xertst( izmin.le.izn .and. izmax.le.izn,
     .                       'Inconsistent charge input!')
      call xertst( imi.gt.0, 'Inconsistent mass input!')
      if(izn .eq. 0) then
        string=' n'
      else if(izn .eq. 1) then
        if(imi .eq. 1) string=' H'
        if(imi .eq. 2) string=' D'
        if(imi .eq. 3) string=' T'
      else
        if(imi.ne.nuclear_mass(izn)) then
          string=' '
          len=1
          if(decorate) string(len+1:len+2)='^{'
          len=lnblnk(string)
          if(imi .lt. 10) then
            write(string(len+1:len+1),'(i1)') imi
          elseif (imi .lt. 100) then
            write(string(len+1:len+2),'(i2)') imi
          else
            write(string(len+1:len+3),'(i3)') imi
          endif
          len=lnblnk(string)
          if(decorate) string(len+1:len+1)='}'
          len=lnblnk(string)
          write(string(len+1:len+lnblnk(elements(izn))),'(a)')
     .     trim(elements(izn))
        else
          string=' ' // trim(elements(izn))
        endif
      endif
      len=lnblnk(string)
      if(decorate) string(len+1:len+2)='^{'
      len=lnblnk(string)
      if(izmin .gt. 0) string(len+1:len+1)='+'
      if(izmin .lt. 0) string(len+1:len+1)='-'
      len=lnblnk(string)
      if(abs(izmin) .lt. 10) then
        write(string(len+1:len+1),'(i1)')
     1   abs(izmin)
      else
        write(string(len+1:len+2),'(i2)')
     1   abs(izmin)
      endif
      len=lnblnk(string)
      if(izmax .ne. izmin) then
        string(len+1:len+1)='_'
        if(izmax .gt. 0) string(len+2:len+2)='+'
        if(izmax .lt. 0) string(len+2:len+2)='-'
        len=lnblnk(string)
        if(abs(izmax) .lt. 10) then
          write(string(len+1:len+1),'(i1)')
     1     abs(izmax)
        else
          write(string(len+1:len+2),'(i2)')
     1     abs(izmax)
        endif
      endif
      len=lnblnk(string)
      if(decorate) string(len+1:len+1)='}'
      call subend()
      return
      end subroutine species
