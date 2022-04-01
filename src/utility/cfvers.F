      subroutine cfvers
      implicit none
#     include <git_version_B25.h>
      character*10 version_in, newversion
      integer nget, nput
      character chf*12,label*7
      logical streql
      external streql

      entry cfverr(nget, newversion)
      inquire (nget,form=chf,err=10)
      if (streql(chf,'formatted')) then
        read(nget,'(a,a)',err=10) label,version_in
      else
        read(nget,err=10) label,version_in
      endif
      if(label.eq.'VERSION') then
	newversion=version_in
        write(*,'(a,a)') 'VERSION set to: ',newversion
        return
      endif
   10 continue
      rewind(nget)
      write(*,'(a,a)') 'VERSION defaulted to : ',newversion
      return

      entry cfverw(nput, newversion)
      inquire (nput,form=chf,err=20)
      if (streql(chf,'formatted')) then
        write(nput,'(a,a,1x,a)') 'VERSION',newversion,git_version_B25
      else
        write(nput) 'VERSION',newversion,' ',git_version_B25
      endif
      write(*,'(4a)') 'VERSION: ',newversion,
     .                ' GIT hash: ',git_version_B25
   20 continue
      return

      end

      function get_B25_hash()
      implicit none
      character*32 get_B25_hash
#     include <git_version_B25.h>

      get_B25_hash = git_version_B25

      return
      end function get_B25_hash

      function get_ADAS_hash()
      implicit none
      character*32 get_ADAS_hash
#     include <git_version_B25.h>

      get_ADAS_hash = git_version_ADAS

      return
      end function get_ADAS_hash

      function get_SOLPS_hash()
      implicit none
      character*32 get_SOLPS_hash
#     include <git_version_B25.h>

      get_SOLPS_hash = git_version_SOLPS

      return
      end function get_SOLPS_hash
