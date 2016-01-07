module logging

  use helper

  ! Loglevels: if ( LEVEL <= LOG... ) then write...
   
  integer, parameter :: LOGFATAL = 1
  integer, parameter :: LOGERROR = 2
  integer, parameter :: LOGWARNING = 3 
  integer, parameter :: LOGKNOWNWARNING = 4
  integer, parameter :: LOGINFO = 5
  integer, parameter :: LOGDEBUG = 6
  integer, parameter :: LOGDEBUGBULK = 7

  integer, parameter :: LOGALL = 1000            ! catchall level

  integer :: loglevel = LOGALL        ! current log level

  logical, parameter :: LOG_WRITELEVEL = .false.  ! write log level of a log message?

contains

  !> Set the given log level
  subroutine setLogLevel( level )
    integer, intent(in) :: level
    
    loglevel = level
  end subroutine setLogLevel

  !> Set the log level given by a specific name
  subroutine setLogLevelByName( levelname )
    character(*), intent(in) :: levelname

    select case ( levelname )
    case ("FATAL", "fatal")
            call setLogLevel(LOGFATAL)
    case ("ERROR", "error")
            call setLogLevel(LOGERROR)
    case ("INFO", "info")
            call setLogLevel(LOGINFO)
    case ("WARNING", "warning")
            call setLogLevel(LOGWARNING)
    case ("KNOWNWARNING", "knownwarning")
            call setLogLevel(LOGKNOWNWARNING)
    case ("DEBUG", "debug")
            call setLogLevel(LOGDEBUG)
    case ("DEBUGBULK", "debugbulk")
            call setLogLevel(LOGDEBUGBULK)
    case default
            stop 'setLogLevelByName: unknown name'
    end select
  end subroutine setLogLevelByName

  ! Returns true if the given log level is encompassed in the current 
  ! log level threshold. Intended usage:
  ! if ( loglvl( LOGDEBUG ) ) write (*,*) ...
  function loglvl( level ) result ( doLog )
    logical :: doLog
    integer, intent(in) :: level

    doLog = ( level <= loglevel ) 
    if ( doLog .and. LOG_WRITELEVEL ) write (*,*) 'loglvl: Loglevel: ', level

  end function loglvl

    
  ! Write log message

  subroutine logmsg( level, msg, ia, ib )
    integer, intent(in) :: level
    character(*), intent(in) :: msg
    integer, intent(in), optional :: ia, ib

    if ( loglvl( level ) ) then
            if ( present( ia ) ) then
                    if ( present( ib ) ) then
                            write (*,'(a,i0,i0)') msg, ia, ib
                    else
                            write (*,'(a,i0)') msg, ia
                    end if
            else 
                    write (*,'(a)') msg
            end if
    end if

  end subroutine logmsg

end module logging
