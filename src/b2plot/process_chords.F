c
      subroutine process_chords(filename)

      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      character*(*) filename
      character*80 dataname,tmpstring
      character*256 extralabel
      logical open_file,vector
#ifndef DEF_MAXCHORDS
#define DEF_MAXCHORDS 10000
#endif
      integer maxchords
      parameter (maxchords=DEF_MAXCHORDS)
      integer i,is,inpts,npoints
      real (kind=R8) ::
     *   xvalues(maxchords),yvalues(maxchords),
     1   xmin,xmax,ymin,ymax,x1,x2,y1,y2,z1,z2
      integer lnblnk
      real (kind=R8) ::
     * chord_integral
      external lnblnk, chord_integral, open_file

      do is=0,ndim(3,iwrk)-1
        if(.not.open_file(99,filename,
     1   'data.local/chords|data/chords')) return
        read(99,*,err=30,end=30) dataname,npoints
        vector=npoints .lt. 0
        if(vector) 
     >   dataname=dataname(1:lnblnk(dataname)) // ' (vector)'
        npoints=abs(npoints)
        inpts=0
   10   read(99,*,err=30,end=40) x1,y1,z1,x2,y2,z2,xvalues(inpts+1)
        inpts=inpts+1
        yvalues(inpts)=chord_integral(nx,ny,work(-1,-1,is,iwrk),
     1   crx,cry,bb,npoints,x1,y1,z1,x2,y2,z2,vector)
        if(mod(inpts,1000) .eq. 0) 
     >   write(*,*) 'Processed ',inpts,' chords'
        if(inpts .lt. maxchords) goto 10
        stop 'increase array size in process_chords'
   30   write(*,*) 'Error in reading ',filename(1:lnblnk(filename))
        return
   40   continue
        close(99)

        if(average .gt. 1) then
          call do_average(xvalues,inpts)
          call do_average(yvalues,inpts)
          inpts=inpts/average
        endif

        xmin=xvalues(1)
        xmax=xvalues(1)
        ymin=yvalues(1)
        ymax=yvalues(1)
        do i=2,inpts
          xmin=min(xmin,xvalues(i))
          xmax=max(xmax,xvalues(i))
          ymin=min(ymin,yvalues(i))
          ymax=max(ymax,yvalues(i))
        enddo
        if(values(11).ne.init_value) xmin=values(11)
        if(values(12).ne.init_value) xmax=values(12)
        if(xmin.eq.xmax) then
          xmin=min(xmin*0.9_R8,xmin*1.1_R8)
          xmax=max(xmax*0.9_R8,xmax*1.1_R8)
          if(xmin.eq.xmax) then
            xmin=-1.0_R8
            xmax=+1.0_R8
          end if
        end if
        if(values(14) .ne. values(15)) then
          ymin=values(14)
          ymax=values(15)
        endif
        if(ymin.eq.ymax) then
          ymin=min(ymin*0.9_R8,ymin*1.1_R8)
          ymax=max(ymax*0.9_R8,ymax*1.1_R8)
          if(ymin.eq.ymax) then
            ymin=-1.0_R8
            ymax=+1.0_R8
          end if
        end if
        tmpstring=textwk(is,iwrk)
        call init_page(
     1   graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk))) // '.m' //
     2   tmpstring(1:lnblnk(tmpstring)) // ' '//
     2   dataname(1:lnblnk(dataname)),
     1   xmin,xmax,ymin,ymax,linlog,.false.,.false.,.false.)
        write(*,*) 'xmin,xmax,ymin,ymax ',xmin,xmax,ymin,ymax
        if(writing_to_file) then
          write(nwtf,'(''# '',a)')                                      !srv 21.09.11
     1     graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk))) // '.m' //
     2     tmpstring(1:lnblnk(tmpstring)) // ' '//
     2     dataname(1:lnblnk(dataname))
          call get_extralabel(extralabel)
          if(len_trim(extralabel).ne.0) then
            write(nwtf,'(''# '',a)') trim(extralabel)                   !srv 21.09.11
          endif
          do i=1,inpts
            write(nwtf,*) xvalues(i),yvalues(i)                         !srv 21.09.11
          enddo
        endif
        do i=1,inpts-1
          call lined_d(xvalues(i),yvalues(i),xvalues(i+1),yvalues(i+1))
        enddo
        call end_page
        close(99)
      enddo
      return
      end subroutine process_chords
c
      subroutine process_camera(filename)

      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
#ifdef CDF_MOVIE
#     include <netcdf.inc>
#endif

      character*(*) filename
      character*80 dataname,tmpstring
#ifndef DEF_MAXCHORDS
#define DEF_MAXCHORDS 10000
#endif
      integer maxchords
      parameter (maxchords=DEF_MAXCHORDS)
      real (kind=R8) ::
     *   zvalues(maxchords)
      integer izvalues(maxchords)
      integer maxcompile
#ifndef HITS_PER_CHORD
#define HITS_PER_CHORD 50
#endif
      parameter (maxcompile=maxchords*HITS_PER_CHORD)
      integer indchd(maxchords+1),indix(maxcompile),indiy(maxcompile)
      real(kind=R8) ::
     * indscl(maxcompile)
#ifdef CDF_MOVIE
      character*5 rw
      integer ncid,imap(maxvdims),iret
      real (kind=R8) ::
     *   htime
      integer itstep,ntstep,dtstep
#endif
      integer inpts,is,i,lnblnk,nxcam,nycam
      external lnblnk
      real (kind=R8) ::
     *   xmin,xmax,ymin,ymax

      if(use_compile) then
        call compile_chords(maxchords,maxcompile,
     &   nx,ny,crx,cry,bb,
     1   inpts,indchd,indix,indiy,indscl,filename,dataname,
     2   nxcam,xmin,xmax,nycam,ymin,ymax)
      endif

#ifdef CDF_MOVIE
      if(moviename .eq. ' ') then
#endif
        do is=0,ndim(3,iwrk)-1
          call calculate_zvalues(maxchords,maxcompile,
     1     inpts,indchd,indix,indiy,indscl,
     1     work(-1,-1,is,iwrk),
     1     zvalues,filename,dataname,nxcam,xmin,xmax,nycam,ymin,ymax)
          tmpstring=textwk(is,iwrk)
          call init_page(
     1     graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk))) //
     2     tmpstring(1:lnblnk(tmpstring)) // ' '//
     2     dataname(1:lnblnk(dataname)),
     1     xmin,xmax,ymin,ymax,1,.true.,.true.,.true.)
          write(*,*) 'xmin,xmax,ymin,ymax ',xmin,xmax,ymin,ymax
          if(writing_to_file) then
            write(nwtf,'(''# '',a)')                   !srv 21.09.11
     1       graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk))) //
     2       tmpstring(1:lnblnk(tmpstring)) // ' '//
     2       dataname(1:lnblnk(dataname))
            do i=1,inpts
              write(nwtf,*) zvalues(i)                   !srv 21.09.11
            enddo
          endif
          call plot_gca(zvalues,izvalues,nxcam,nycam,
     1     xmin,ymin,xmax,ymax,values(9),values(10),
     2     linlog)
          call end_page
        enddo
#ifdef CDF_MOVIE
      else
        rw='read'
        ncid=ncopn(moviename,ncnowrit,iret)
        imap(1)=1
        call rwcdf(rw,ncid,'ntstep',imap,ntstep,iret)
        do itstep=1,ntstep,movie_skip
          call rwcdf_settime('time',itstep)
          call rwcdf(rw,ncid,'time',imap,htime,iret)
          imap(1)=1
          imap(2)=nx+2
          call rwcdf(rw,ncid,'ahal',imap,work(-1,-1,0,iwrk),iret)
          call calculate_zvalues(maxchords,maxcompile,
     1     inpts,indchd,indix,indiy,indscl,
     1     work(-1,-1,0,iwrk),
     1     zvalues,filename,dataname,nxcam,xmin,xmax,nycam,ymin,ymax)
          write(tmpstring,'(f10.3,a)')
     1     htime*1.0e6_R8,
     2     ' us'
          call init_page(
     1     graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk))) //
     2     tmpstring(1:lnblnk(tmpstring)),
     1     xmin,xmax,ymin,ymax,1,.true.,.true.,.true.)
          write(*,*) 'xmin,xmax,ymin,ymax ',xmin,xmax,ymin,ymax
          if(writing_to_file) then
            write(nwtf,'(''# '',a)')                   !srv 21.09.11
     1       graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk))) //
     2       tmpstring(1:lnblnk(tmpstring)) // ' '//
     2       dataname(1:lnblnk(dataname))
            do i=1,inpts
              write(nwtf,*) zvalues(i)                   !srv 21.09.11
            enddo
          endif
          call plot_gca(zvalues,izvalues,nxcam,nycam,
     1     xmin,ymin,xmax,ymax,values(9),values(10),
     2     linlog)
          call end_page
        enddo
        call ncclos(ncid, iret)
      endif
#endif
      return
      end subroutine process_camera
c
      subroutine compile_chords(maxchords,maxcompile,
     & nx,ny,crx,cry,bb,
     1 inpts,indchd,indix,indiy,indscl,filename,dataname,
     2 nxcam,xmin,xmax,nycam,ymin,ymax)

      use b2mod_types
      implicit none

      integer maxchords,maxcompile,nx,ny,inpts
      integer indchd(maxchords+1),indix(maxcompile),indiy(maxcompile)
      real (kind=R8) ::
     * indscl(maxcompile),
     1 crx(-1:nx,-1:ny,0:3),
     2 cry(-1:nx,-1:ny,0:3),
     3 bb(-1:nx,-1:ny,0:3)
      character*(*) filename,dataname
      integer nxcam,nycam
      real (kind=R8) ::
     *   xmin,ymin,xmax,ymax

      logical open_file,vector
      integer ind,indold,npoints,lnblnk
      external lnblnk
      real (kind=R8) ::
     *   x1,x2,y1,y2,z1,z2

      if(.not.open_file(99,filename,
     1 'data.local/chords|data/chords')) return
      read(99,*,err=30,end=30) dataname,npoints,
     1 nxcam,xmin,xmax,nycam,ymin,ymax
      vector=npoints .lt. 0
      if(vector)
     > dataname=dataname(1:lnblnk(dataname)) // ' (vector)'
      npoints=abs(npoints)
      inpts=0
      ind=1
   10 read(99,*,err=30,end=40) x1,y1,z1,x2,y2,z2
      inpts=inpts+1
      call chord_compile(npoints,nx,ny,crx,cry,bb,
     1 x1,y1,z1,x2,y2,z2,vector,
     1 indchd(inpts),indix(ind),indiy(ind),indscl(ind))
      indold=ind
      ind=ind+indchd(inpts)
      if(ind .gt. maxcompile)
     1 stop 'increase maxcompile in process_camera'
      indchd(inpts)=indold
      if(mod(inpts,1000) .eq. 0) 
     > write(*,*) 'Processed ',inpts,' chords'
      if(inpts .lt. maxchords) goto 10
      stop 'increase array size in process_chords'
   30 write(*,*) 'Error in reading ',filename(1:lnblnk(filename))
      return
   40 continue
      indchd(inpts+1)=ind
      write(*,*) ind-1,' elements of ',maxcompile,' used'
      close(99)
      if(inpts .ne. nxcam*nycam) then
        write(*,*) 'Problem with chords, ',nxcam*nycam,
     1   ' expected, ',inpts,' read in'
        return
      endif
      return
      end subroutine compile_chords
c
      subroutine calculate_zvalues(maxchords,maxcompile,
     1 inpts,indchd,indix,indiy,indscl,workarr,
     2 zvalues,filename,dataname,
     3 nxcam,xmin,xmax,nycam,ymin,ymax)

      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot
      implicit none

      integer maxchords,maxcompile,inpts
      integer indchd(maxchords+1),indix(maxcompile),indiy(maxcompile)
      real (kind=R8) ::
     * indscl(maxcompile)
      real (kind=R8) ::
     *   zvalues(maxchords)
      real (kind=R8) ::
     *   workarr(-1:nx,-1:ny)
      character*(*) filename,dataname
      integer nxcam,nycam
      real (kind=R8) ::
     *   xmin,ymin,xmax,ymax

      logical open_file,vector
      integer npoints,lnblnk,i,j
      real (kind=R8) ::
     *   x1,x2,y1,y2,z1,z2,chord_integral
      external lnblnk,chord_integral

      if(use_compile) then
        do i=1,inpts
          zvalues(i)=0.0_R8
          do j=indchd(i),indchd(i+1)-1
            zvalues(i)=zvalues(i)+
     1       workarr(indix(j),indiy(j))*indscl(j)
          enddo
        enddo
      else
        if(.not.open_file(99,filename,
     1   'data.local/chords|data/chords')) return
        read(99,*,err=30,end=30) dataname,npoints,
     1   nxcam,xmin,xmax,nycam,ymin,ymax
        vector=npoints .lt. 0
        if(vector)
     1   dataname=dataname(1:lnblnk(dataname)) // ' (vector)'
        npoints=abs(npoints)
        inpts=0
   10   read(99,*,err=30,end=40) x1,y1,z1,x2,y2,z2
        inpts=inpts+1
        zvalues(inpts)=chord_integral(nx,ny,workarr(-1,-1),
     1   crx,cry,bb,npoints,x1,y1,z1,x2,y2,z2,vector)
        if(mod(inpts,1000) .eq. 0)
     >   write(*,*) 'Processed ',inpts,' chords'
        if(inpts .lt. maxchords) goto 10
        stop 'increase array size in process_chords'
   30   write(*,*) 'Error in reading ',filename(1:lnblnk(filename))
        return
   40   continue
        close(99)
        if(inpts .ne. nxcam*nycam) then
          write(*,*) 'Problem with chords, ',nxcam*nycam,
     1     ' expected, ',inpts,' read in'
          return
        endif
      endif
      return
      end subroutine calculate_zvalues
