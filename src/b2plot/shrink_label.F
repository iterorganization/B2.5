
      subroutine shrink_label(label)
      implicit none
      integer lng, i, idot, ie
      character*(*) label

      lng=len_trim(label)
C     write(*,*) ' In: ',trim(label)
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
      lng=len_trim(label)
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
      if(label.eq.' '.or.label.eq.'.') then
        label(1:1)='0'
      elseif(label(lng:lng).eq.'.') then
        label(lng:lng)=' '
      endif
C     write(*,*) 'Out: ',trim(label)
      return
      end subroutine shrink_label
