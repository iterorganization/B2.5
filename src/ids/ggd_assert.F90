module ggd_assert

  !> Assertions for Fortran 90, written by H.-J. Klingshirn

  !> @author 

  use b2mod_types
  use string 

  implicit none

  private

  public assertSetStopMode, assertSetMsgPrefix, assertStopOnFailed, assert, assertEqual, assertReset

  logical, save :: defaultDoStop = .true.
  integer, save :: failCount = 0
  integer, save :: msgPrefixLen = 0
  character(256), save :: msgPrefix

contains

  !> Private subroutine that handles a failed assertion. Prints out the
  !> fail message (msgPrefix + failmsg, depending on what is given), and 
  !> either stops of does assertion bookkeeping. 

  !> @param failmsg The fail message to print
  !> @param doStop Controls whether to stop execution. If doStop .true., 
  !> the program is stopped. If .false., only the fail message is printed
  !> and bookkeeping is done for delayed stopping (@see assertStopOnFailed)

  subroutine assertFail( failmsg, doStop )
    character(*), intent(in), optional :: failmsg
    logical, intent(in), optional :: doStop

    ! internal
    logical :: lDoStop

    lDoStop = defaultDoStop
    if ( present( doStop ) ) lDoStop = doStop

    if ( present( failmsg ) ) then
            if ( msgPrefixLen > 0 ) then
                    write (*,*) 'assert: '//msgPrefix(1:msgPrefixLen)//' '//failmsg
            else
                    write (*,*) 'assert: '//failmsg       
            end if
    else
            if ( msgPrefixLen > 0 ) then
                    write (*,*) 'assert: '//msgPrefix(1:msgPrefixLen)
            end if
    end if

    if ( lDoStop ) then
       stop "assert: an assertion failed, stopping immediately."
    else
       failCount = failCount + 1
       write(*,*) "assert: an assertion failed, continuing."
    end if
    
  end subroutine assertFail


  !> Set the default assertion stop behaviour. 

  !> @param doStop .true. means immediate
  !> stop on fail, .false. means don't stop on fail (for use with 
  !> delayed stopping, @see assertStopOnFailed). The default value is .true.,
  !> it is set if doStop is omitted.
  subroutine assertSetStopMode( doStop )
    logical, intent(in), optional :: doStop

    if ( present( doStop ) ) then
            defaultDoStop = doStop
    else
            ! if no argument given, just set default again
            defaultDoStop = .true.           
    end if
    
  end subroutine assertSetStopMode


  !> Set a prefix for the assertion fail message.
  !> @param prefix The prefix string for the fail messages. Optional. If 
  !> omitted, clear the prefix string (i.e. no prefix is used anymore)
  subroutine assertSetMsgPrefix( prefix )
    character(len=*), intent(in), optional :: prefix

    if ( present( prefix ) ) then
            msgPrefixLen = min( len(prefix), len( msgPrefix ) )
            msgPrefix(1:msgPrefixLen) = prefix(1:msgPrefixLen)
    else
            msgPrefixLen = 0
    end if
    
  end subroutine assertSetMsgPrefix


  !> Reset the assertion module to its default state.
  !>
  !> Reset stop mode and message prefix to their default values.
  !> Also resets the fail counter, i.e. any previous failed assertions are forgotton.
  subroutine assertReset()
    call assertSetStopMode()       
    call assertSetMsgPrefix()
    failCount = 0
  end subroutine assertReset
  

  !> Stop if a previously called assertion failed (delayed stop). 

  !> This is useful if the
  !> default stop mode is set to continue on failed assertions (which
  !> can be enabled by call assertSetStopMode( .false. ). 

  !> @param failmsg Message to print if assertion(s) failed.
  !> @param doStop Controls stop behaviour, @see assert. Overrides the default
  !> set with @see assertSetStopMode

  subroutine assertStopOnFailed( failmsg, doStop )
    character(*), intent(in), optional :: failmsg
    logical, intent(in), optional :: doStop

    ! internal
    logical :: lDoStop 
    
    lDoStop = .true. ! default for this routine is to stop, regardless of current module default    
    if ( present( doStop ) ) lDoStop = doStop ! ...but can be overridden again with optional argument

    if ( present( failmsg ) ) then
            ! ...use given message
            call assert( failCount == 0, failmsg, lDoStop )
    else
            ! ...or substitute generic message
            call assert( failCount == 0, Int2Str( failCount )//' assertions failed', lDoStop )
    end if
    failCount = 0

  end subroutine assertStopOnFailed
  

  !> A generic assertion, tests a given logical expression. 

  !> If it evaluates to .false.,
  !> print the fail message and possibly stop execution.
  !> 
  !> @param test The logical expression to test.
  !> @param failmsg The message to print on fail. If omitted, a generic message is printed. Can be modified with a prefix (see assertSetMsgPrefix)
  !> @param doStop Controls whether to stop execution. If doStop .true., 
  !> the program is stopped. If .false., only the fail message is printed
  !> and bookkeeping is done for delayed stopping (see assertStopOnFailed).
  !> If given, overrides the default behaviour  set by assertSetStopMode.

  !> @see assertStopOnFailed
  !> @see assertSetStopMode
  !> @see assertSetMsgPrefix

  !> @author H.-J. Klingshirn 
  !> @version 1.0 

  subroutine assert( test, failmsg, doStop )
    logical, intent(in) :: test
    character(*), intent(in), optional :: failmsg
    logical, intent(in), optional :: doStop

    if ( .not. test ) then
       call assertFail( failmsg, doStop )
    end if

  end subroutine assert


  !> Test double precision floating point numbers for equality.

  !> @param x1, x2 The values to test
  !> @param failmsg Same as for assert
  !> @param doStop Same as for assert
  !> @see assert

  subroutine assertEqual( x1, x2, failmsg, doStop )
    real(R8), intent(in) :: x1, x2
    character(*), intent(in), optional :: failmsg
    logical, intent(in), optional :: doStop

    ! x1 is reference value
    ! x2 is actual value

    if ( .not. ( ( ( x1 + 2 * spacing( real( x1 ) ) ) >= x2 ) &
           & .and. ( x1 - 2 * spacing( real( x1 ) ) ) <= x2 ) ) then

       call assertFail( failmsg, doStop )

    end if
    
  end subroutine assertEqual


end module ggd_assert

!!!Local Variables:
!!! mode: f90
!!! End:
