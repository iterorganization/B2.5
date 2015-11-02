      subroutine cfvers
      implicit none
#     include <git_version.h>
      character*10 version_in, newversion
      integer nget, nput
      character chf*12,label*7
      logical, save :: unread
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
        write(*,*) 'VERSION set to: ',newversion
        return
      endif
   10 continue
      rewind(nget)
      write(*,*) 'VERSION defaulted to : ',newversion
      return

      entry cfverw(nput, newversion)
      inquire (nput,form=chf,err=20)
      if (streql(chf,'formatted')) then
        write(nput,'(a,a,1x,a)') 'VERSION',newversion,git_version
      else
        write(nput) 'VERSION',newversion,' ',git_version
      endif
      write(*,'(4a)') 'VERSION: ',newversion,' GIT hash: ',git_version
   20 continue
      return

      end
