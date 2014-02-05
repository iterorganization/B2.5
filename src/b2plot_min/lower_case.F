c
      subroutine lower_case(sval)                   !srv 12.11.12
      use b2mod_types
      character*(*) sval
      integer i
      character*1 c
      do i=1,len(sval)
        c=sval(i:i)
        if(c .ge. 'A' .and. c .le. 'Z') then
          sval(i:i)=char(ichar(c)-ichar('A')+ichar('a'))
        endif
      enddo
      return
      end
