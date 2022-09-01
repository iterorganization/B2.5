      subroutine strip_spaces(string)
      implicit none
      integer :: i,j,n
      character*1 tab
      character*(*) string

      tab=char(9)
      j=0
      n=len(string)
      do i=1,n
        if(string(i:i).ne.' '.and.string(i:i).ne.tab) then
          j=j+1
          if(i.ne.j) string(j:j)=string(i:i)
        endif
      enddo
      if(j.lt.n) string(j+1:n)=' '
      return
      end subroutine strip_spaces
