! Compare output files from two different B2.5 runs.
! It can be used to compare b2fxxx files
!
! Compile: ifort -g -O2 check_b2_output.F90 -o check_b2_output
!
! Usage: ./check_b2_output original_run/b2mn.exe.dir/b2fstate run2/b2mn.exe.dir/b2fstate 

module b2mod_types_local
  implicit none
  INTEGER, PARAMETER :: R8 = SELECTED_REAL_KIND (14)
  INTEGER, PARAMETER :: R4 = SELECTED_REAL_KIND (6)
end module b2mod_types_local


module check_module_local
  implicit none
  private
  public :: check_variable

  double precision, parameter :: tolerance = 1.0d-14

  ! compare two arrays, and reports if there is any difference
  interface check_variable
    module procedure check_variable_c0, check_variable_i1, check_variable_r1
  end interface

 contains 

 function check_variable_c0(val, orig, name) result(n_error)
  character(len=*), intent(in) :: val, orig
  character(len=*), intent(in) :: name
  integer :: n_error 
  if (val.ne.orig) then
    write (*,*) 'Error in ',name,', values: ', val, orig
    n_error = 1
  else
    n_error = 0
  end if
  end function 

  function check_variable_i1(val, orig, name) result(n_error)
  integer, dimension(:), intent(in) :: val, orig
  character(len=*), intent(in) :: name
  integer :: n_error
  integer :: error, error_max
  integer :: i1
  
  n_error = 0
  error_max = 0
  if (ubound(val,1).ne.ubound(orig,1)) then
    write (*,*) 'Error: array boundary differ ', name, 1, ubound(val,1), ubound(orig,1)
  endif

  do i1 = 1, ubound(val, 1)
    error = abs(val(i1) - orig(i1))
    if (error > tolerance) then
      if (n_error < 5) then
        write (*,*) 'Error in ',name,  error
        write (*,*) '  values: ', val(i1), orig(i1)
        write (*,*) '  at indices: ',i1
      end if
      n_error = n_error + 1
      if (error > error_max) error_max = error
    end if
  end do
  if (n_error > 0) then
    if (n_error.eq.1) then
      write (*,*) 'There was 1 error in array ', name
    else
      write (*,*) 'There were ', n_error, 'errors in array ', name
    endif
    write (*,*) 'Max error is ', error_max, ' in array ', name
  end if
  end function

  function check_variable_r1(val, orig, name) result(n_error)
  double precision, dimension(:), intent(in) :: val, orig
  character(len=*), intent(in) :: name
  integer :: n_error

  double precision :: error, error_max, avg_abs_diff, avg_rel_error, avg_abs_val
  double precision :: error_max_abs
  integer :: i1, ub, max_err_idx
  
  n_error = 0
  error_max = 0
  error_max_abs = 0
  ub = ubound(val,1)
  if (ubound(val,1).ne.ubound(orig,1)) then
      write (*,*) 'Error: array boundary differ ', name, 1, ubound(val,1), ubound(orig,1)
      if ( ubound(val,1) > ubound(orig,1)) then
        ub = ubound(orig,1)
      endif
  endif

  avg_abs_val = 0
  avg_abs_diff = 0
  avg_rel_error = 0

  do i1 = 1, ub
    if (orig(i1)==0 .or. val(i1)==0) then
      error = abs(val(i1) - orig(i1))
    else
      error = abs((val(i1) - orig(i1)) / orig(i1))
    end if
    avg_abs_val = avg_abs_val + abs(orig(i1))
    if (error > tolerance) then
      if (n_error < 5) then
        write (*,*) 'Error in ',name,  error
        write (*,*) '  values: ', val(i1), orig(i1)
        write (*,*) '  at index: ',i1
      end if
      n_error = n_error + 1
      avg_abs_diff = avg_abs_diff + abs(val(i1) - orig(i1))
      avg_rel_error = avg_rel_error + error
      if (error > error_max) then
         error_max = error
         error_max_abs =  abs(val(i1) - orig(i1))
         max_err_idx = i1
      endif
    end if
  end do
  if (n_error > 0) then
    if (n_error.eq.1) then
      write (*,*) 'There was 1 error in array ', trim(name)
    else
      write (*,*) 'There were ', n_error, 'errors in array ', trim(name)
    end if
    write (*,'(a13,a10, E14.5E3,a8,i8)') 'Max error in ', trim(name), error_max, ' at idx ' , max_err_idx
    write (*,*) 'values', orig(max_err_idx), val(max_err_idx)
    write (*,'(a22,a10, E14.5E3)') 'Max absolute error in ', trim(name), error_max_abs
    write (*,'(a28,a10,E14.5E3)') 'Average relative error in    ',trim(name), avg_rel_error / n_error
    write (*,'(a28,a10,E14.5E3)') 'Average absolute diff in     ',trim(name), avg_abs_diff / n_error
    avg_abs_val = avg_abs_val / ub
    write (*,'(a29,a10,E14.5E3)') 'Average array value (abs) in ',trim(name), avg_abs_val 
    write (*,'(a29,a10,E14.5E3)') 'Max error / avg array value  ',trim(name), error_max_abs / avg_abs_val 
    write (*,'(a30,F5.1,a1)')   'Number of errors / array size ', n_error * 100.0 / ub, '%'
    write (*,*) ' '
  end if
end function
end module

module b2_file_io
! Code pieces extracted from B2
  implicit none
  private
  public read_unknown_type !< read in an array when we do not know which type to expect
 
  contains
 
  subroutine read_record_description(nget, chf, idcod, idtyp, n, id, ierr)
      integer, intent(in) :: nget
      character, intent(out) :: chf*12, idcod*8, idtyp*8, id*32
      integer, intent(out) :: n
      integer, optional, intent(out) :: ierr
      integer :: ierror
      inquire (nget,form=chf)
      if (chf=='formatted' .or. chf=='FORMATTED') then
        chf = 'FORMATTED'
        read (nget,'(2a8,i12,4x,a32)',iostat=ierror) idcod, idtyp, n, id
      else
        read (nget,iostat=ierror) idcod, idtyp, n, id
      endif
      if (present(ierr)) ierr = ierror
  end subroutine

  subroutine check_record_description(idcod, id0, id1, idtyp0, idtyp1, n0, n1, strict)
    character(len=*), intent(in) :: idcod, idtyp0, idtyp1, id0, id1
    integer, intent(in) :: n0, n1
    logical, optional :: strict
    if ('*cf:' /= idcod) then
      write (*,*)  'cfrure--not a header line, expected *cf:, found '//idcod
      stop
    endif
    if (idtyp0 /= idtyp1) then
      write(*,*)   'cfrure--wrong type: expected '//idtyp0//', found '//idtyp1
      stop
    endif
    if (present(strict)) then
      if (strict) then
        if (n0 /= n1) then
          write(*,*) 'cfrure--dimension mismatch at '//id0
          write(*,*) 'Expected size: ', n0, 'Found: ', n1
          stop
        endif
        if (id0 /= id1) then
          write (*,*) 'cfrure--wrong label: expected '//id0//', found '//id1
          stop
        endif
      endif
    endif
  end subroutine

  subroutine read_unknown_type(nget, n1, refun, infun, chfun, id1, idtyp)
!     The file stores real, integer or character data.
!     We check which type of data is next, and read it into one 
!     of the output buffers (refun, infun or chfun)
!     The variable name (id1) and type (idtyp) are also returned,
!     together with the size of the array that is read (n1)
    use b2mod_types_local   !IGNORE
    implicit none
      integer, intent(in) :: nget  !< file descriptor
      integer, intent(out) :: n1 !< array size
      !> Buffers to store real, integer or character data
      real(kind=r8), dimension(:), allocatable, intent(out) :: refun
      integer, dimension(:), allocatable, intent(out) :: infun
      character(len=:), allocatable, intent(out) :: chfun
      character(len=32), intent(inout) :: id1 !< variable name
      character(len=8), intent(out) :: idtyp  !< variable type
      integer :: i, ierr, n
      character chf*12, idcod*8
!     ------------------------------------------------------------------
      call read_record_description(nget, chf, idcod, idtyp, n1, id1, ierr)
      if (ierr < 0) then
         n1 = 0
         id1 = 'end-of-file'
         return
      endif
      call check_record_description(idcod, id1, id1, idtyp, idtyp, n1, n1)
      select case (idtyp)
        case ('real')
          allocate(refun(0:n1-1))
          if (0.lt.n1) then
            if (chf=='FORMATTED') then
              read (nget,*) (refun(i),i=0,n1-1)
            else
              read (nget) (refun(i),i=0,n1-1)
            endif
          endif  
        case ('int')
          allocate(infun(0:n1-1))
          if (0.lt.n1) then
            if (chf=='FORMATTED') then
              read (nget,*) (infun(i),i=0,n1-1)
            else
              read (nget) (infun(i),i=0,n1-1)
            endif
          endif
        case ('char')
          n = n1
          allocate(character(n) :: chfun)
          if (0.lt.n1) then
            if (chf=='FORMATTED') then
              read (nget,'(1x,a)') chfun(1:n1)
            else
              read (nget) chfun(1:n1)
            endif
          endif
     end select
  end subroutine
end module


program test_b2output
  use check_module_local   !IGNORE
  use b2mod_types_local    !IGNORE
  use b2_file_io           !IGNORE
  implicit none

  character(len=200) :: input1, input2
  integer :: u1, u2, idx
  real (kind=r8), dimension(:), allocatable :: r1, r2
  integer, dimension(:), allocatable :: i1, i2
  character(len=:), allocatable :: ch1, ch2
  integer :: size_n1, size_n2, n_errors
  character(len=32) :: vname1, vname2
  character(len=8) :: idtyp1, idtyp2

  call get_filenames(input1, input2)
  u1 = open_file(input1)
  u2 = open_file(input2)

  n_errors = 0 
  call read_unknown_type (u1, size_n1, r1, i1, ch1, vname1, idtyp1)
  call read_unknown_type (u2, size_n2, r2, i2, ch2, vname2, idtyp2)
  
  ! Read variables (arrays) one by one from the files, and compare them.
  do while(vname1 /= 'end-of-file')
    if (vname1 /= vname2 .or. (idtyp1 /= idtyp2)) then
      write(*,*) 'Error, variable name or type differ'
      write(*,*) trim(vname1), ' ', idtyp1
      write(*,*) trim(vname2), ' ', idtyp2
    else 
        
        select case(idtyp1)
          case ('real')
            write(*,'(a15,a10,a7,i8,a6,a8)') 'Checking array ', vname1, ' size: ', size_n1, ' type ', idtyp1
            n_errors = n_errors + check_variable(r1, r2, vname1)
          case ('int')
             write(*,'(a15,a10,a7,i8,a6,a8)') 'Checking array ', vname1, ' size: ', size_n1, ' type ', idtyp1
            n_errors = n_errors + check_variable(i1, i2, vname1)
          case ('char')
             write(*,'(a16,a10,a6,i8)') 'Skipping string ', vname1, 'size: ', size_n1 
             ! might contain time when the simulation was started
            !n_errors = n_errors + check_variable(ch1, ch2, vname1)
          case default
            write(*,*) 'unknown type', idtyp1, ' for variable ', vname1
        end select
    endif
    call read_unknown_type (u1, size_n1, r1, i1, ch1, vname1, idtyp1)
    call read_unknown_type (u2, size_n2, r2, i2, ch2, vname2, idtyp2)
  enddo

  write(*,*)
  idx = index(input1, '/', back=.true.) + 1
  if (n_errors == 0) then
    write (*,*) 'There were no errors in ', input1(idx:len_trim(input1))
  else if (n_errors == 1) then
    write (*,*) 'There was 1 error in ', input1(idx:len_trim(input1))
  else
    write (*,*) 'There were', n_errors, 'errors in ', input1(idx:len_trim(input1))
  endif

  close(u1)
  close(u2)
  
 contains

 subroutine get_filenames(input1, input2)
  ! Get the filenames from the command line arguments
  character(len=200), intent(out) :: input1, input2
  logical :: ex
  call get_command_argument(1, input1)
  call get_command_argument(2, input2)
  if (LEN_TRIM(input1) == 0 .or. LEN_TRIM(input2)==0) then
    write (*,*) 'Error, two input files shoud be given as arguments!'
    stop
  endif

  inquire(file=input1, exist=ex)
  if (.not. ex) then
    write (*,*) 'Error: input file not found:', input1
    stop
  endif
  inquire(file=input2, exist=ex)
  if (.not. ex) then
    write (*,*) 'Error: input file not found:', input2
    stop
  endif

  write (*,*) 'File1: ', trim(input1)
  write (*,*) 'File2: ', trim(input2)
 end subroutine
 
 function open_file(filename) result(my_unit)
! Open the file, read the header, and return the file unit specifier.
! B2.5 output files can be either binary or text files. We open first as
! formatted file, if it does not make sense, then reopen as unformatted.
   character(*), intent(in) :: filename
   logical:: use_formatted_io
   integer :: my_unit
   integer, dimension(:), allocatable :: idum
   character(len=:), allocatable ::  lblmn
   character(len=32) :: id
   integer :: n 
   character(len=10) :: version_in
   character(len=7) :: label
   integer :: ierr
   logical :: b2fstate
   
   open(newunit=my_unit,file=trim(filename), status='old', action='read', form='FORMATTED', iostat=ierr)
   ! Read the header
   read(my_unit,'(a,a)') label ,version_in
   if (label/='VERSION') then
     if (label == '*cf:') then
       ! No version header, but seems to be a text file
       label='unknown'
       version_in = 'version'
       rewind(my_unit)
     else
       ! If it does not match, then it should be a binary file
       close(my_unit)
       ! reopen in UNFORMATTED mode
       open(newunit=my_unit,file=trim(filename), status='old', action='read', form='UNFORMATTED', iostat=ierr)
       read(my_unit) label ,version_in
       if (label/='VERSION') then
         if (label == '*cf:') then
           ! No version header, but seems to be a valid file
           label='unknown'
           version_in = 'version'
           rewind(my_unit)
         else
           write(*,*) 'Error, file should start with VERSION or *cf tag'
           stop
         endif
       endif
     endif
   endif
   if (ierr == 0) write(*,*) trim(filename),' opened successfully'
   write(*,*) label, ' ', version_in
 end function
end program

!!!Local Variables:
!!! mode: f90
!!! End:
