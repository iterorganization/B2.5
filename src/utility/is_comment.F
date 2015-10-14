c
      logical function is_comment(line)
      character*(*) line
      logical streql
      external streql
      is_comment = streql(line,' ') .or.
     .             line(1:1).eq.'*' .or.
     .             line(1:1).eq.'!' .or.
     .             line(1:1).eq.'#'
      return
      end
