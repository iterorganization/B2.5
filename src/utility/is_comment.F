c
      logical function is_comment(line)
      character*(*), intent(in) :: line
      logical streql
      external streql
      is_comment = streql(line,' ') .or.
     .             line(1:1).eq.'*' .or.
     .             line(1:1).eq.'!' .or.
     .             line(1:1).eq.'#'
      return
      end function is_comment

      logical function is_valid(cval)
      character*(*), intent(in) :: cval
      character*80 :: cvbl
      logical streql, isadigit, isaletter, isplusorminus
      external streql, isadigit, isaletter, isplusorminus

      cvbl = cval
      call strip_spaces(cvbl)
      is_valid = .false.
      if (.not.streql(cvbl,' '))
     > is_valid = isadigit(cvbl(1:1)) .or. isaletter(cvbl(1:1)) .or.
     &  (isplusorminus(cvbl(1:1)).and.isadigit(cvbl(2:2)))
      return
      end function is_valid

      logical function end_comment(line,ind)
      character*(*), intent(in) :: line
      integer, intent(out) :: ind
      logical streql
      external streql
      integer i, j, k, l

      i = index(line,'''')
      j = 0
      if (i.gt.0) j = index(line(i+1:),'''')
      k = 0
      if (j.gt.0) k = index(line(i+j+1:),'''')
      l = 0
      if (k.gt.0) l = index(line(i+j+k+1:),'''')
      ind = 0
      if (l.gt.0) ind = verify(line(i+j+k+l+1:),' ')
      end_comment = ind.gt.0
      ind = ind + i + j + k + l

      return
      end function end_comment

      logical function only_comments(line)
      character*(*), intent(in) :: line
      integer i

      only_comments = .true.
      do i = 1, len(line)
        only_comments = only_comments .and.
     &   (line(i:i).eq.' '.or.line(i:i).eq.'*'.or.
     &    line(i:i).eq.'!'.or.line(i:i).eq.'#')
      end do

      return
      end function only_comments
