      module b2mod_solpstop

      contains

      function solpstop()
#ifdef NAGF95
      use f90_unix           ! IGNORE
#endif
      implicit none
#ifndef NO_GETENV
#ifndef USE_PXFGETENV
      integer lenval, ierror
      intrinsic get_environment_variable
#endif
#endif
      character*256 solpstop, save_solpstop
      logical first
      save first,save_solpstop
      data first /.true./

      if(first) then
#ifdef NO_GETENV
        save_solpstop='/afs/ipp-garching.mpg.de/u/dpc/solps/solps-iter'
#else
#ifdef USE_PXFGETENV
        CALL PXFGETENV ('SOLPSTOP', 0, save_solpstop, lenval, ierror)
#else
        call get_environment_variable('SOLPSTOP',
     .   status=ierror, length=lenval)
        if (ierror.eq.0) then
          call get_environment_variable('SOLPSTOP',value=save_solpstop)
        elseif (ierror.eq.1) then
          save_solpstop=
     .     '/afs/ipp-garching.mpg.de/u/dpc/solps/solps-iter'
        endif
#endif
#endif
        write(*,*) 'SOLPSTOP = : ',trim(save_solpstop)
        first=.false.
      endif
      solpstop=save_solpstop
      return

      end function solpstop

      end module b2mod_solpstop
