      subroutine strip_spaces(string)
      Implicit None
      integer :: i,j,n
      character*(*) string
      j=0
      n=len(string)
      do i=1,n
        if(string(i:i).ne.' ') then
          j=j+1
          if(i.ne.j) string(j:j)=string(i:i)
        endif
      enddo
      if(j.lt.n) string(j+1:n)=' '
      return
      end
