c
      subroutine vesselon(onflag)
      use b2mod_types
#include "graphics_vars.h"
      logical onflag
      logical exist
      character*256 filename

      if(onflag) then
        vessel=.not.vessel
        write (*,*) 'vessel state changed'
        filename='mesh.extra'
        call find_file(filename,exist)
        if (.not.exist.and.vessel) then
          write(*,*) 'No mesh.extra file found !'
          vessel=.false.
          write(*,*) 'Unable to use vessel option.'
        end if
      else
        vessel=.false.
        write (*,*) 'vessel off'
      endif
      return
      end subroutine vesselon
