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
      logical first
c
      save first, solpstop
      data first /.true./
c
      if(first) then
#ifdef NO_GETENV
        solpstop=' '
#else
#ifdef USE_PXFGETENV
        CALL PXFGETENV ('SOLPSTOP', 0, solpstop, lenval, ierror)
#else
        call getenv ('SOLPSTOP', solpstop)
#endif
#endif
        if(solpstop.eq.' ') then
          write(*,*) 'open_file: no SOLPSTOP set '
        else
          write(*,*) 'open_file: SOLPSTOP = : ',
     1	   solpstop(1:lnblnk(solpstop))
        endif
        first=.false.
      endif
      open_file=.true.
      close(lu)
      write(*,*) 'Trying ',filename(1:lnblnk(filename))
      open(lu,file=filename(1:lnblnk(filename)),
     1 status='old',err=10,recl=1000010)
      return
 10   continue
      close(lu)
      filenamecat='../'//filename(1:lnblnk(filename))
      write(*,*) 'Trying ',filenamecat(1:lnblnk(filenamecat))
      open(lu,file=filenamecat(1:lnblnk(filenamecat)),
     1 status='old',err=12,recl=1000010)             !srv 12.05.08 {
      return
 12   continue                           
      close(lu)
      filenamecat='/data/adas/adf15/'//filename(1:lnblnk(filename))
      write(*,*) 'Trying ',filenamecat(1:lnblnk(filenamecat))
      open(lu,file=filenamecat(1:lnblnk(filenamecat)),
     1 status='old',err=13,recl=1000010)             !srv 12.05.08 } !srv 15.07.09 {
      return
 13   continue                           
      close(lu)
      filenamecat='/data/'//filename(1:lnblnk(filename))
      write(*,*) 'Trying ',filenamecat(1:lnblnk(filenamecat))
      open(lu,file=filenamecat(1:lnblnk(filenamecat)),
     1 status='old',err=15,recl=1000010)                             !srv 15.07.09 }
      return
 15   continue
      copyarea=area
 17   iend=index(copyarea,'|')
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
      write(*,*) 'Trying ',filenamecat(1:lnblnk(filenamecat))
      open(lu,file=filenamecat(1:lnblnk(filenamecat)),
     3 status='old',err=20)
      return
 20   continue
      if(iend.gt.0) goto 17
      write(*,*) 'Error in opening ',filename(1:lnblnk(filename))
      open_file=.false.
      return
c
      end

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
      logical first
c
      save first, solpstop
      data first /.true./
c
      if(first) then
#ifdef NO_GETENV
        solpstop=' '
#else
#ifdef USE_PXFGETENV
        CALL PXFGETENV ('SOLPSTOP', 0, solpstop, lenval, ierror)
#else
        call getenv ('SOLPSTOP', solpstop)
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
 10   continue
      close(lu)
      filenamecat='../'//filename(1:lnblnk(filename))
      open(lu,file=filenamecat(1:lnblnk(filenamecat)),
     1 status='old',err=15,recl=1000010)
      return
 15   continue
      copyarea=area
 17   iend=index(copyarea,'|')
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
 20   continue
      if(iend.gt.0) goto 17
      silent_open=.false.
      write(*,*) 'Failed opening ', filename(1:lnblnk(filename))
      write(*,*) 'SOLPSTOP = ', solpstop(1:lnblnk(solpstop))
      write(*,*) 'area = ', area(1:lnblnk(area))
      filenamecat=''
      return
c
      end
