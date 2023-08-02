!nh   20.06.2023
!     Subroutine copied from EIRENE and slightly altered. 
!     Necessary to read in reaction data in b2ar.dat.    
      subroutine read_token (inchar, divide, outchar, itok, ier,
     .  lreal)
 
      implicit none
 
      character(len=*), intent(in) :: inchar, divide
      character(len=*), intent(out) :: outchar
      integer, intent(out) :: itok, ier
      logical, intent(in) :: lreal
 
      integer :: len_in, len_out, ianf, iend, len_tok, iend2
 
! determine length of character strings
      len_in = len(inchar)
      len_out = len(outchar)
 
! initialize outchar
      outchar = repeat(' ',len_out)
 
! determine position of first non-divide character
      ianf = verify(inchar,divide)
 
      if (ianf > 0) then
! a non-divide character has been found
        itok = ianf-1 + scan(inchar(ianf:),divide)
        iend = itok-1
        if (lreal) then
          if ((len(divide) == 1) .and. (divide(1:1) == ' ')) then
            if (index('eEdD',inchar(iend:iend)) > 0) then
              iend2 = scan(inchar(itok+1:),divide)
              iend = iend + iend2
              itok = itok + iend2
            end if
          end if
        end if
        len_tok = iend-ianf+1
 
        ier = 0
        if (len_tok > len_out) then
! string too long for output variable
          ier = 1
          outchar(1:len_out) = inchar(ianf:ianf+len_out-1)
        else
! string fits into output variable
          outchar(1:len_tok) = inchar(ianf:iend)
        end if
      else
! string contains only divide-characters
        ier = 0
        itok = 0
        outchar(1:1) = '0'
      end if
 
      return
      end subroutine read_token
 
 
 
