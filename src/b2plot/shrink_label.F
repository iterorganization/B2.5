
      subroutine shrink_label(label)
      integer lnblnk, lng, i, idot, ie
      external lnblnk
      character*(*) label
      lng=lnblnk(label)
C     write(*,*) ' In: ',label(1:lnblnk(label))
   10 continue
      do i=1,lng
        if(label(i:i).eq.' ') then
c remove spaces
          label(i:lng)=label(min(i+1,lng):lng)//' '
          lng=lng-1
          goto 10
        elseif(label(i:i+1).eq.'0E') then
c x0Ey -> xEy
          label(i:lng)=label(min(i+1,lng):lng)//' '
          lng=lng-1
          goto 10
        elseif(label(i:i+4).eq.'.E+00') then
c x.E+00 -> x
          label(i:lng)='     '
          goto 20
        elseif(label(i:i+3).eq.'.E+0') then
c x.E+0y -> xEy
          label(i:lng)='E'//label(min(i+4,lng):lng)//'    '
          goto 20
        elseif(label(i:i+2).eq.'.E+') then
c x.E+y -> xEy
          label(i:lng)='E'//label(min(i+3,lng):lng)//'    '
          goto 20
        elseif(label(i:i+3).eq.'.E-0') then
c x.E-0y -> xE-y
          label(i:lng)='E-'//label(min(i+4,lng):lng)//'    '
          goto 20
        elseif(label(i:i+2).eq.'.E-') then
c x.E-y -> xE-y
          label(i:lng)='E-'//label(min(i+3,lng):lng)//'    '
          goto 20
        elseif(label(i:i+3).eq.'E+00') then
c xE+00 -> x
          label(i:lng)='     '
          goto 20
        elseif(label(i:i+2).eq.'E+0') then
c xE+0y -> xEy
          label(i:lng)='E'//label(min(i+3,lng):lng)//'    '
          goto 20
        elseif(label(i:i+1).eq.'E+') then
c xE+y -> xEy
          label(i:lng)='E'//label(min(i+2,lng):lng)//'    '
          goto 20
        elseif(label(i:i+2).eq.'E-0') then
c xE-0y -> xE-y
          label(i:lng)='E-'//label(min(i+3,lng):lng)//'    '
          goto 20
        endif
      enddo
   20 continue
      lng=lnblnk(label)
      idot=0
      ie=0
      do i=1,lng
        if(label(i:i).eq.'.') idot=i
        if(label(i:i).eq.'E') ie=i
      enddo
      if(idot.ne.0.and.ie.eq.0) then
        do while(label(lng:lng).eq.'0')
          label(lng:lng)=' '
          lng=lng-1
        enddo
      endif
      if(label(lng:lng).eq.'.') label(lng:lng)=' '
      if(label.eq.' ') label(1:1)='0'
C      write(*,*) 'Out: ',label(1:lnblnk(label))
      return
      end
