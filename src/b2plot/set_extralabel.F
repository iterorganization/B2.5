
C ========== Parse extralabel for customisation =====================
      subroutine customise_label(locallabel,customlabel)
C ===================================================================
      use b2mod_types
      use b2mod_b2plot_debug
#ifdef NAGF95
      use f90_unix           ! IGNORE
      use f90_unix_proc      ! IGNORE
#endif
#include <graphics_vars.h>
#include <COMLAB_local.F>
      character*(*) locallabel
      character*256 tmp,customlabel
      integer elength,count,i,j,ijpos,icpos,jcpos,istart,iend
      integer lnblnk
      logical streql
      external lnblnk, subini, subend, xertst, xerrab, streql
      character var*256,envvar*256,systemarg*256
      character*1 backslash
#ifdef USE_PXFGETENV
      integer lenval, ierror
#else
#ifdef NAGFOR
      integer lenval, ierror
#endif
#endif
*!!   (We avoid to declare, as external or intrinsic, the subroutine
*     getenv; it is external on some systems (SunOS, Unicos, ...) and
*     intrinsic on others (AIX, HP-UX, ...).  If getenv is declared
*     external on AIX or HP-UX it returns some non-printing garbage,
*     without error message (AIX/XLF Version 2.2, HP-UX Release 9.0))
*     ------------------------------------------------------------------
      call subini ('customise_label')
      backslash='\\'
      icpos=1
      jcpos=1
      elength=lnblnk(locallabel)
      if(streql(locallabel,' ')) goto 1000
#ifndef PLPLOT
      call pcsetc('FC - function code character','@')
#endif
   10 istart=icpos
      do i=istart,elength
        if(locallabel(icpos:icpos).eq.backslash) then
          if(locallabel(icpos+1:icpos+1).eq.'R') goto 12
          if(locallabel(icpos+1:icpos+1).eq.'G') goto 14
          if(locallabel(icpos+1:icpos+1).eq.'d') goto 16
          if(locallabel(icpos+1:icpos+1).eq.'g') goto 20
          if(locallabel(icpos+1:icpos+1).eq.'r') goto 30
          if(locallabel(icpos+1:icpos+1).eq.'n') goto 40
          if(locallabel(icpos+1:icpos+1).eq.'h') goto 50
          if(locallabel(icpos+1:icpos+1).eq.'t') goto 60
        endif
        if(locallabel(icpos:icpos+1).eq.'^{') goto 80
        if(locallabel(icpos:icpos+1).eq.'_{') goto 90
        if(locallabel(icpos:icpos+1).eq.'${') goto 100
        if(locallabel(icpos:icpos).eq.'`') goto 110
        tmp(jcpos:jcpos)=locallabel(icpos:icpos)
        icpos=icpos+1
        jcpos=jcpos+1
      enddo
      goto 1000
   12 continue ! Switch to the "Roman" part of the PWRITX database
      tmp(jcpos:jcpos+2)='@R@'
      icpos=icpos+2
      jcpos=jcpos+3
      goto 10
   14 continue ! Switch to the "Greek" part of the PWRITX database
      tmp(jcpos:jcpos+2)='@G@'
      icpos=icpos+2
      jcpos=jcpos+3
      goto 10
   16 continue ! default characters
      tmp(jcpos:jcpos+3)='@F0@'
      icpos=icpos+2
      jcpos=jcpos+4
      goto 10
   20 continue ! found greek letters
      tmp(jcpos:jcpos+4)='@F33@'
      icpos=icpos+2
      jcpos=jcpos+5
      goto 10
   30 continue ! found roman letters
      tmp(jcpos:jcpos+4)='@F25@'
      icpos=icpos+2
      jcpos=jcpos+5
      goto 10
   40 continue ! switch to normal level
      tmp(jcpos:jcpos+2)='@N@'
      icpos=icpos+2
      jcpos=jcpos+3
      goto 10
   50 continue ! resize font 130%
      tmp(jcpos:jcpos+5)='@Z150@'
      icpos=icpos+2
      jcpos=jcpos+6
      goto 10
   60 continue ! resize font 70%
      tmp(jcpos:jcpos+4)='@Z50@'
      icpos=icpos+2
      jcpos=jcpos+5
      goto 10
   80 continue ! found superscript
      tmp(jcpos:jcpos+2)='@S@'
      icpos=icpos+1
      jcpos=jcpos+3
      if(locallabel(icpos:icpos).eq.'{') then
        icpos=icpos+1
        ijpos=icpos
        do j=ijpos,elength
          if(locallabel(j:j).eq.'}') then
            tmp(jcpos:jcpos+2)='@N@'
            icpos=icpos+1
            jcpos=jcpos+3
            goto 10
          else
            tmp(jcpos:jcpos)=locallabel(icpos:icpos)
            icpos=icpos+1
            jcpos=jcpos+1
          endif
        enddo
        call xertst (j.lt.elength,
     .   'You did not close superscript braces')
      endif
      goto 10
   90 continue ! found subscript
      tmp(jcpos:jcpos+2)='@B@'
      icpos=icpos+1
      jcpos=jcpos+3
      if(locallabel(icpos:icpos).eq.'{') then
        icpos=icpos+1
        ijpos=icpos
        do j=ijpos,elength
          if(locallabel(j:j).eq.'}') then
            tmp(jcpos:jcpos+2)='@N@'
            icpos=icpos+1
            jcpos=jcpos+3
            goto 10
          else
            tmp(jcpos:jcpos)=locallabel(icpos:icpos)
            icpos=icpos+1
            jcpos=jcpos+1
          endif
        enddo
        call xertst(j.lt.elength, 'You did not close subscript braces')
      endif
      goto 10
  100 continue ! found locallabel
      icpos=icpos+2
      iend=icpos
      do while (iend.lt.elength.and.locallabel(iend:iend).ne.'}')
        iend=iend+1
      enddo
      call xertst(iend.lt.elength.or.locallabel(iend:iend).eq.'}',
     . 'You did not close environment braces')
      var=locallabel(icpos:iend-1)
      if(var.eq.'time'.or.var.eq.'TIME') then
        if((runtime.lt.1e-6.and.runtime.ne.0).or.(runtime.ge.100)) then
          write(envvar,'(1pg11.4)') runtime
        else
          write(envvar,'(f9.6)') runtime
        endif
      elseif(var.eq.'label'.or.var.eq.'LABEL') then
        envvar=label
      elseif(var.eq.'sonnetfile'.or.var.eq.'SONNETFILE') then
        envvar=sonnetfile
      else
#ifdef NO_GETENV
        envvar=' '
#else
#ifdef NAGFOR
        call get_environment_variable(var, status=ierror, length=lenval)
        if (ierror.eq.0) then
          call get_environment_variable(var,value=envvar)
        elseif (ierror.eq.1) then
          envvar=' '
        endif
#else
#ifdef USE_PXFGETENV
        CALL PXFGETENV (var, 0, envvar, lenval, ierror)
#else
        call getenv (var, envvar)
#endif
#endif
#endif
        if(envvar.eq.' ') then
          write(*,*) 'Environment variable ',var(1:lnblnk(var)),
     1       ' does not exist'
        endif
      endif
      tmp(jcpos:jcpos+lnblnk(envvar)-1)=envvar(1:lnblnk(envvar))
      jcpos=jcpos+lnblnk(envvar)
      icpos=iend+1
      goto 10
  110 continue ! `
      icpos=icpos+1
      iend=icpos
      do while (iend.lt.elength.and.locallabel(iend:iend).ne.'`')
        iend=iend+1
      enddo
      call xertst(iend.lt.elength.or.locallabel(iend:iend).eq.'`',
     1 'You did not close `')
      systemarg='echo `'//locallabel(icpos:iend-1)//'` > /tmp/tmp.xxx'
      call system(systemarg)
      open(99,file='/tmp/tmp.xxx')
      read(99,'(a)') envvar
      close(99)
      call system('/bin/rm /tmp/tmp.xxx')
      tmp(jcpos:jcpos+lnblnk(envvar)-1)=envvar(1:lnblnk(envvar))
      jcpos=jcpos+lnblnk(envvar)
      icpos=iend+1
      goto 10
 1000 continue
      if(b2plot_debug.ge.10) write (*,*) 'jcpos =',jcpos
      if(jcpos.lt.256) then
        do i=jcpos,256
          tmp(i:i)=' '
        enddo
      else
        write(*,*) 'jcpos = ',jcpos
        write(*,*) 'tmp = ',tmp
!        call xerrab ('extralabel too long')
        jcpos=256
      endif
      count=0
      do i=1,jcpos
        if(tmp(i:i).eq.'@'.and.tmp(i+1:i+1).eq.'@') then
          count=count+2
          do j=i,jcpos-2
            tmp(j:j)=tmp(j+2:j+2)
          enddo
        endif
      enddo
      jcpos=jcpos-count
      do i=jcpos,256
        tmp(i:i)=' '
      enddo
      if(b2plot_debug.ge.4) then
        if (streql(tmp,' ')) then
          write (*,*) 'parsed empty extralabel'
        else
          write (*,*) 'parsed extralabel to:   ',trim(tmp)
        end if
      end if
      customlabel=tmp
      call subend ()
      return
      end subroutine customise_label
c
      subroutine get_extralabel(return_extralabel)
C ===================================================================
      use b2mod_types
      use b2mod_b2plot_debug
#include <graphics_vars.h>
#include <COMLAB_local.F>
      character*(*) return_extralabel
      return_extralabel=extralabel
      return
      end subroutine get_extralabel
c
      subroutine set_wallbarlabel(locallabel)
C ===================================================================
      implicit none
#include <graphics_vars.h>
      character*(*) locallabel
      integer lnblnk
      external lnblnk
      wallbarlabel=locallabel(1:lnblnk(locallabel))
      return
      end subroutine set_wallbarlabel
c
      subroutine set_targetbarlabel(locallabel)
C ===================================================================
      implicit none
#include <graphics_vars.h>
      character*(*) locallabel
      integer lnblnk
      external lnblnk
      targetbarlabel=locallabel(1:lnblnk(locallabel))
      return
      end subroutine set_targetbarlabel
