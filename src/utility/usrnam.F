      character*16 function usrnam ()
      use b2mod_types
#ifdef NAGF95
      use f90_unix           ! IGNORE
#endif
      implicit none
*     ------------------------------------------------------------------
*     USRNAM returns the name or other appropriate identification of
*     the submitter of the present job.
*!!   This routine may contain system dependencies.
*     This version employs the system routine getenv to acquire the
*     value of an appropriate environment variable.
*     ------------------------------------------------------------------
      character usrnm*16
#ifndef NO_GETENV
      integer lenval, ierror
#ifndef USE_PXFGETENV
      intrinsic get_environment_variable
#endif
#endif
      save usrnm
      data usrnm/' '/
*!!   (We avoid to declare, as external or intrinsic, the subroutine
*     getenv; it is external on some systems (SunOS, Unicos, ...) and
*     intrinsic on others (AIX, HP-UX, ...).  If getenv is declared
*     external on AIX or HP-UX it returns some non-printing garbage,
*     without error message (AIX/XLF Version 2.2, HP-UX Release 9.0))
*     ------------------------------------------------------------------
      if (usrnm.eq.' ') then
*      (this is the first call to usrnam)
#ifndef NO_GETENV
#ifdef USE_PXFGETENV
       CALL PXFGETENV ('LOGNAME', 0, usrnm, lenval, ierror)
#else
       call get_environment_variable('LOGNAME',
     .   status=ierror, length=lenval)
       if (ierror.eq.0)
     >   call get_environment_variable('LOGNAME', value=usrnm)
#endif
       if (usrnm.eq.' ') then
#ifdef USE_PXFGETENV
        CALL PXFGETENV ('USER', 0, usrnm, lenval, ierror)
#else
        call get_environment_variable('USER',
     .    status=ierror, length=lenval)
        if (ierror.eq.0)
     >    call get_environment_variable('USER', value=usrnm)
#endif
        if (usrnm.eq.' ') then
#ifdef USE_PXFGETENV
         CALL PXFGETENV ('user', 0, usrnm, lenval, ierror)
#else
         call get_environment_variable('user',
     .     status=ierror, length=lenval)
         if (ierror.eq.0)
     >     call get_environment_variable('user', value=usrnm)
#endif
#endif
         if (usrnm.eq.' ') then
           usrnm = 'anonymous'
         endif
        endif
       endif
      endif
      usrnam = usrnm
      return
*     ------------------------------------------------------------------
      end function usrnam
