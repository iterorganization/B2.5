      logical function open_file(lu,filename,area)
      use b2mod_types
c
#ifdef NAGF95
      use f90_unix           ! IGNORE
#endif
      implicit none
c
      integer lu
#ifndef NO_GETENV
      integer lenval, ierror
#ifndef USE_PXFGETENV
      intrinsic get_environment_variable
#endif
#endif
*!!   (We avoid to declare, as external or intrinsic, the subroutine
*     getenv; it is external on some systems (SunOS, Unicos, ...) and
*     intrinsic on others (AIX, HP-UX, ...).  If getenv is declared
*     external on AIX or HP-UX it returns some non-printing garbage,
*     without error message (AIX/XLF Version 2.2, HP-UX Release 9.0))
*     ------------------------------------------------------------------
      character*(*) filename,area
c
      character*256 solpstop,filenamecat,localarea,copyarea
      integer iend
      logical streql
      external streql
      logical first
c
      save first, solpstop
      data first /.true./
c
      if(first) then
        solpstop=' '
#ifndef NO_GETENV
#ifdef USE_PXFGETENV
        CALL PXFGETENV ('SOLPSTOP', 0, solpstop, lenval, ierror)
#else
        call get_environment_variable('SOLPSTOP',
     .   status=ierror, length=lenval)
        if (ierror.eq.0) call get_environment_variable('SOLPSTOP',
     .   value=solpstop)
#endif
#endif
        if(streql(solpstop,' ')) then
          write(*,'(a)') 'open_file: no SOLPSTOP set '
        else
          write(*,'(a,a)') 'open_file: SOLPSTOP = : ',
     1     trim(solpstop)
        endif
        first=.false.
      endif
      open_file=.true.
      close(lu)
      write(*,'(a,a)') 'Trying ',trim(filename)
      open(lu,file=trim(filename),status='old',err=10,recl=1000010)
      return
   10 continue
      close(lu)
      filenamecat='../'//trim(filename)
      write(*,'(a,a)') 'Trying ',trim(filenamecat)
      open(lu,file=trim(filenamecat),status='old',err=12,recl=1000010)   !srv 12.05.08 {
      return
   12 continue
      close(lu)
      filenamecat=trim(solpstop)//'/'//trim(filename)
      write(*,'(a,a)') 'Trying ',trim(filenamecat)
      open(lu,file=trim(filenamecat),status='old',err=13)
      return
   13 continue
      close(lu)
      copyarea=area
      if (.not.streql(copyarea,' ')) then
   17   iend=index(copyarea,'|')
        if(iend.eq.0) then
          localarea=copyarea
        else
          localarea=copyarea(1:iend-1)
          copyarea(1:256-iend)=copyarea(iend+1:256)
          copyarea(256-iend+1:256)=repeat(' ',iend)
        endif
        filenamecat=
     1   trim(solpstop)//'/'//trim(localarea)//'/'//trim(filename)
        write(*,'(a,a)') 'Trying ',trim(filenamecat)
        open(lu,file=trim(filenamecat),status='old',err=20)
        return
   20   continue
        close(lu)
        if(iend.gt.0) goto 17
      endif
      write(*,'(a,a)') 'Error in opening ',trim(filename)
      open_file=.false.
      return
c
      end function open_file

      logical function silent_open(lu,filename,area,filenamecat)
      use b2mod_types
c
#ifdef NAGF95
      use f90_unix           ! IGNORE
#endif
      implicit none
c
      integer lu
#ifndef NO_GETENV
      integer lenval, ierror
#ifndef USE_PXFGETENV
      intrinsic get_environment_variable
#endif
#endif
*!!   (We avoid to declare, as external or intrinsic, the subroutine
*     getenv; it is external on some systems (SunOS, Unicos, ...) and
*     intrinsic on others (AIX, HP-UX, ...).  If getenv is declared
*     external on AIX or HP-UX it returns some non-printing garbage,
*     without error message (AIX/XLF Version 2.2, HP-UX Release 9.0))
*     ------------------------------------------------------------------
      character*(*) filename,area
c
      character*256 solpstop,filenamecat,localarea,copyarea
      integer iend
      logical first, streql
      external streql
c
      save first, solpstop
      data first /.true./
c
      if(first) then
        solpstop=' '
#ifndef NO_GETENV
#ifdef USE_PXFGETENV
        CALL PXFGETENV ('SOLPSTOP', 0, solpstop, lenval, ierror)
#else
        call get_environment_variable('SOLPSTOP',
     .   status=ierror, length=lenval)
        if (ierror.eq.0) call get_environment_variable('SOLPSTOP',
     .   value=solpstop)
#endif
#endif
        first=.false.
      endif
      silent_open=.true.
      close(lu)
      open(lu,file=trim(filename),status='old',err=10,recl=1000010)
      filenamecat=trim(filename)
      return
   10 continue
      close(lu)
      filenamecat='../'//trim(filename)
      open(lu,file=trim(filenamecat),status='old',err=12,recl=1000010)
      return
   12 continue
      close(lu)
      filenamecat='../baserun/'//trim(filename)
      open(lu,file=trim(filenamecat),status='old',err=14,recl=1000010)
      return
   14 continue
      close(lu)
      filenamecat='../../baserun/'//trim(filename)
      open(lu,file=trim(filenamecat),status='old',err=16,recl=1000010)
      return
   16 continue
      copyarea=area
   17 iend=index(copyarea,'|')
      if(iend.eq.0) then
         localarea=copyarea
      else
         localarea=copyarea(1:iend-1)
         copyarea(1:256)=copyarea(iend+1:256)
      endif
      close(lu)
      if (.not.streql(localarea,' ')) then
        filenamecat=
     1   trim(solpstop)//'/'//trim(localarea)//'/'//trim(filename)
      else
        filenamecat=
     1   trim(solpstop)//'/'//trim(filename)
      endif
      open(lu,file=trim(filenamecat),status='old',err=20)
      return
   20 continue
      if(iend.gt.0) goto 17
      silent_open=.false.
      write(*,'(a,a)') 'Failed opening ', trim(filename)
      write(*,'(a,a)') 'SOLPSTOP = ', trim(solpstop)
      if(.not.streql(area,'')) write(*,'(a,a)') 'area = ', trim(area)
      filenamecat=''
      return
c
      end function silent_open

      subroutine find_file(filename,exist)
      implicit none

! this routine checks whether the file exists in the present directory,
! one level up, or in baserun (assumed parallel to the directory one
! level up) and returns a new filename if found, and exist set to true
! if the file is found

      character(len=*) filename
      logical exist

      inquire(file=filename, exist=exist)
      if(exist) return

      inquire(file='../'//filename, exist=exist)
      if(exist) then
         filename='../'//trim(filename)
         return
      endif

      inquire(file='../baserun/'//filename, exist=exist)
      if(exist) then
         filename='../baserun/'//trim(filename)
         return
      endif

      inquire(file='../../baserun/'//filename, exist=exist)
      if(exist) then
         filename='../../baserun/'//trim(filename)
         return
      endif

      return
      end subroutine find_file
