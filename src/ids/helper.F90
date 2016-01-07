module helper

  implicit none

contains  

  function int2str(intval) result(string)
    integer , intent(in) :: intval
    character(len=2-max(sign(1,intval),0)+max( &
         min(abs(intval)/10**1,1)*1, &
         min(abs(intval)/10**2,1)*2, &
         min(abs(intval)/10**3,1)*3, &
         min(abs(intval)/10**4,1)*4, &
         min(abs(intval)/10**5,1)*5, &
         min(abs(intval)/10**6,1)*6, &
         min(abs(intval)/10**7,1)*7, &
         min(abs(intval)/10**8,1)*8, &
         min(abs(intval)/10**9,1)*9) ) :: string
    integer :: absn,j,k,is
    absn = abs(intval)
    if ( absn == intval ) then
            is = 1
    else
            is = 2
            string(1:1) = "-"
    end if
    do j=len(string),is,-1
            k = modulo(absn,10)+1
            string(j:j) = "0123456789"(k:k)
            absn = absn / 10
    end do
    return
  end function int2str

  !> Write a real to a string with 15 digits
  !> +1.23456789012345E000
  function real2str(realval) result(string)
    character(*), parameter :: SAMPLE = '+1.23456789012345E000'
    !                                     0.10000000000000E+03
    !                                    -0.10000000000000E+03
    double precision, intent(in) :: realval
    character(len(SAMPLE)) :: string

    write( string, '(es21.14)' ) realval
  end function real2str

end module helper
