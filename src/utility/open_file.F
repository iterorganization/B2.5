      logical function open_file(lu,filename,area)
      use b2mod_types
c
#ifdef NAGF95
      use f90_unix           ! IGNORE
#endif
      implicit none
c
      integer lu
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
      character*(*) filename,area
c
      character*256 solpstop,filenamecat,localarea,copyarea
      integer lnblnk, iend
      logical streql
      external lnblnk, streql
      logical first
c
      save first, solpstop
      data first /.true./
c
      if(first) then
#ifdef NO_GETENV
        solpstop=' '
#else
#ifdef NAGFOR
        call get_environment_variable('SOLPSTOP',
     .   status=ierror, length=lenval)
        if (ierror.eq.0) then
          call get_environment_variable('SOLPSTOP',value=solpstop)
        elseif (ierror.eq.1) then
          solpstop=' '
        endif
#else
#ifdef USE_PXFGETENV
        CALL PXFGETENV ('SOLPSTOP', 0, solpstop, lenval, ierror)
#else
        call getenv ('SOLPSTOP', solpstop)
#endif
#endif
#endif
        if(solpstop.eq.' ') then
          write(*,'(a)') 'open_file: no SOLPSTOP set '
        else
          write(*,'(a,a)') 'open_file: SOLPSTOP = : ',
     1     solpstop(1:lnblnk(solpstop))
        endif
        first=.false.
      endif
      open_file=.true.
      close(lu)
      write(*,'(a,a)') 'Trying ',filename(1:lnblnk(filename))
      open(lu,file=filename(1:lnblnk(filename)),
     1 status='old',err=10,recl=1000010)
      return
   10 continue
      close(lu)
      filenamecat='../'//filename(1:lnblnk(filename))
      write(*,'(a,a)') 'Trying ',filenamecat(1:lnblnk(filenamecat))
      open(lu,file=filenamecat(1:lnblnk(filenamecat)),
     2 status='old',err=12,recl=1000010)             !srv 12.05.08 {
      return
   12 continue
      close(lu)
      filenamecat=solpstop(1:lnblnk(solpstop))//
     1 '/'//filename(1:lnblnk(filename))
      write(*,'(a,a)') 'Trying ',filenamecat(1:lnblnk(filenamecat))
      open(lu,file=filenamecat(1:lnblnk(filenamecat)),
     3 status='old',err=13)
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
        filenamecat=solpstop(1:lnblnk(solpstop))//
     1   '/'//
     2   localarea(1:lnblnk(localarea))//
     3   '/'//
     4   filename(1:lnblnk(filename))
        write(*,'(a,a)') 'Trying ',filenamecat(1:lnblnk(filenamecat))
        open(lu,file=filenamecat(1:lnblnk(filenamecat)),
     4   status='old',err=20)
        return
   20   continue
        close(lu)
        if(iend.gt.0) goto 17
      endif
      write(*,'(a,a)') 'Error in opening ',filename(1:lnblnk(filename))
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
      character*(*) filename,area
c
      character*256 solpstop,filenamecat,localarea,copyarea
      integer lnblnk, iend
      logical first, streql
      external streql
c
      save first, solpstop
      data first /.true./
c
      if(first) then
#ifdef NO_GETENV
        solpstop=' '
#else
#ifdef NAGFOR
        call get_environment_variable('SOLPSTOP',
     .   status=ierror, length=lenval)
        if (ierror.eq.0) then
          call get_environment_variable('SOLPSTOP',value=solpstop)
        elseif (ierror.eq.1) then
          solpstop=' '
        endif
#else
#ifdef USE_PXFGETENV
        CALL PXFGETENV ('SOLPSTOP', 0, solpstop, lenval, ierror)
#else
        call getenv ('SOLPSTOP', solpstop)
#endif
#endif
#endif
        first=.false.
      endif
      silent_open=.true.
      close(lu)
      open(lu,file=filename(1:lnblnk(filename)),
     1 status='old',err=10,recl=1000010)
      filenamecat=filename(1:lnblnk(filename))
      return
   10 continue
      close(lu)
      filenamecat='../'//filename(1:lnblnk(filename))
      open(lu,file=filenamecat(1:lnblnk(filenamecat)),
     1 status='old',err=12,recl=1000010)
      return
   12 continue
      close(lu)
      filenamecat='../baserun/'//filename(1:lnblnk(filename))
      open(lu,file=filenamecat(1:lnblnk(filenamecat)),
     1 status='old',err=14,recl=1000010)
      return
   14 continue
      close(lu)
      filenamecat='../../baserun/'//filename(1:lnblnk(filename))
      open(lu,file=filenamecat(1:lnblnk(filenamecat)),
     1 status='old',err=16,recl=1000010)
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
      if (lnblnk(localarea).gt.0) then
        filenamecat=solpstop(1:lnblnk(solpstop))//
     1   '/'//
     2   localarea(1:lnblnk(localarea))//
     3   '/'//
     4   filename(1:lnblnk(filename))
      else
        filenamecat=solpstop(1:lnblnk(solpstop))//
     1   '/'//
     2   filename(1:lnblnk(filename))
      endif
      open(lu,file=filenamecat(1:lnblnk(filenamecat)),
     3 status='old',err=20)
      return
   20 continue
      if(iend.gt.0) goto 17
      silent_open=.false.
      write(*,'(a,a)') 'Failed opening ', filename(1:lnblnk(filename))
      write(*,'(a,a)') 'SOLPSTOP = ', solpstop(1:lnblnk(solpstop))
      if(.not.streql(area,''))
     > write(*,'(a,a)') 'area = ', area(1:lnblnk(area))
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

      integer lnblnk
      character(len=*) filename
      logical exist
      external lnblnk

      inquire(file=filename, exist=exist)
      if(exist) return

      inquire(file='../'//filename, exist=exist)
      if(exist) then
         filename='../'//filename(1:lnblnk(filename))
         return
      endif

      inquire(file='../baserun/'//filename, exist=exist)
      if(exist) then
         filename='../baserun/'//filename(1:lnblnk(filename))
         return
      endif

      inquire(file='../../baserun/'//filename, exist=exist)
      if(exist) then
         filename='../../baserun/'//filename(1:lnblnk(filename))
         return
      endif

      return
      end subroutine find_file

