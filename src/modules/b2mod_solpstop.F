      module b2mod_solpstop

      contains

      function solpstop()
#ifdef NAGF95
      use f90_unix           ! IGNORE
#endif
      implicit none
#ifdef USE_PXFGETENV
      integer lenval, ierror
#endif
      character*256 solpstop, save_solpstop
      logical first
      save first,save_solpstop
      data first /.true./

      if(first) then
#ifdef NO_GETENV
	save_solpstop='/afs/ipp-garching.mpg.de/u/dpc/solps/solps5.0'
#else
#ifdef USE_PXFGETENV
        CALL PXFGETENV ('SOLPSTOP', 0, save_solpstop, lenval, ierror)
#else
        call getenv ('SOLPSTOP', save_solpstop)
#endif
#endif
        write(*,*) 'SOLPSTOP = : ',trim(save_solpstop)
        first=.false.
      endif
      solpstop=save_solpstop
      return

      end function solpstop

      end module b2mod_solpstop
