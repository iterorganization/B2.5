c
      subroutine process_eirene_chords(filename)

      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot
      use b2mod_b2plot_debug
      use b2mod_b2plot_values
#ifdef B25_EIRENE
      use eirmod_ctrig
      implicit none
      character*(*) filename
      character*80 dataname
      character*256 extralabel,tmpstring,titlestring
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
     * chord_tria_integral
      external chord_tria_integral
      external lnblnk, open_file, xerrab, subini, subend

      call subini ('process_eirene_chords')
      do is=0,nedim(2,iewrk)-1
        if(.not.open_file(99,filename,
     1   'data.local/chords|data/chords')) then
          call subend ()
          return
        end if
        read(99,*,err=30,end=30) dataname,npoints
        vector=npoints .lt. 0
        if(vector) 
     >   dataname=dataname(1:lnblnk(dataname)) // ' (vector)'
        npoints=abs(npoints)
        inpts=0
   10   read(99,*,err=30,end=40) x1,y1,z1,x2,y2,z2,xvalues(inpts+1)
        inpts=inpts+1
        yvalues(inpts)=chord_tria_integral(ntrii,nrknot,nx,ny,
     1   eirene_stack(1,is,iewrk),
     2   xtrian,ytrian,necke,ixtri,iytri,crx,cry,bb,npoints,
     3   x1,y1,z1,x2,y2,z2,vector)
        if(mod(inpts,1000) .eq. 0) 
     >   write(*,*) 'Processed ',inpts,' chords'
        if(inpts .lt. maxchords) goto 10
        call xerrab ('increase array size in process_eirene_chords')
   30   write(*,*) 'Error in reading ',filename(1:lnblnk(filename))
        call subend ()
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
        if(values(11).ne.init_value .and.
     &      .not.(values(11).eq.0.0_R8 .and. values(12).eq.0.0_R8)) 
     >   xmin=values(11)
        if(values(12).ne.init_value .and.
     &      .not.(values(11).eq.0.0_R8 .and. values(12).eq.0.0_R8)) 
     >   xmax=values(12)
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
        tmpstring=
     1   eirene_label(iewrk)(1:lnblnk(eirene_label(iewrk))) // '.m ' //
     2   textewk(is,iewrk)(1:lnblnk(textewk(is,iewrk))) // ' ' //
     3   dataname(1:lnblnk(dataname))
        call customise_label(tmpstring,titlestring)
        call init_page(titlestring,
     1   xmin,xmax,ymin,ymax,linlog,.false.,.false.,.false.)
        if(b2plot_debug.ge.4)
     >   write(*,*) 'xmin,xmax,ymin,ymax ',xmin,xmax,ymin,ymax
        if(writing_to_file) then
          write(nwtf,'(''# '',a)') titlestring(1:lnblnk(titlestring))
          call get_extralabel(extralabel)
          if(len_trim(extralabel).ne.0) then
            write(nwtf,'(''# '',a)') trim(extralabel)
          endif
          do i=1,inpts
            write(nwtf,*) xvalues(i),yvalues(i)
          enddo
        endif
        do i=1,inpts-1
          call lined_d(xvalues(i),yvalues(i),xvalues(i+1),yvalues(i+1))
        enddo
        call end_page
        close(99)
      enddo
      call subend ()
#endif
      return
      end subroutine process_eirene_chords
c
      subroutine process_eirene_camera(filename)

      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot
      use b2mod_b2plot_debug
      use b2mod_b2plot_values
#ifdef B25_EIRENE
      use eirmod_ctrig
#endif
      implicit none

      character*(*) filename
#ifdef B25_EIRENE
      character*80 dataname
      character*256 tmpstring,titlestring
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
      integer indchd(maxchords+1),indit(maxcompile)
      real(kind=R8) ::
     * indscl(maxcompile)
      integer inpts,is,i,lnblnk,nxcam,nycam
      external lnblnk
      real (kind=R8) ::
     *   xmin,xmax,ymin,ymax

      if(use_compile) then
        call compile_tria_chords(maxchords,maxcompile,
     &   ntrii,nrknot,nx,ny,
     &   xtrian,ytrian,necke,ixtri,iytri,crx,cry,bb,
     1   inpts,indchd,indit,indscl,filename,dataname,
     2   nxcam,xmin,xmax,nycam,ymin,ymax)
      endif

      do is=0,nedim(2,iewrk)-1
        call calculate_tzvalues(maxchords,maxcompile,
     1   inpts,indchd,indit,indscl,
     1   eirene_stack(1,is,iewrk),
     1   zvalues,filename,dataname,nxcam,xmin,xmax,nycam,ymin,ymax)
        tmpstring=
     1   eirene_label(iewrk)(1:lnblnk(eirene_label(iewrk))) //
     2   textewk(is,iewrk)(1:lnblnk(textewk(is,iewrk))) // ' ' //
     3   dataname(1:lnblnk(dataname))
        call customise_label(tmpstring,titlestring)
        call init_page(titlestring,
     1   xmin,xmax,ymin,ymax,1,.true.,.true.,.true.)
        if(b2plot_debug.ge.4)
     >   write(*,*) 'xmin,xmax,ymin,ymax ',xmin,xmax,ymin,ymax
        if(writing_to_file) then
          write(nwtf,'(''# '',a)') titlestring(1:lnblnk(titlestring))
          do i=1,inpts
            write(nwtf,*) zvalues(i)
          enddo
        endif
        call plot_gca(zvalues,izvalues,nxcam,nycam,
     1   xmin,ymin,xmax,ymax,values(9),values(10),
     2   linlog)
        call end_page
      enddo
#endif
      return
      end subroutine process_eirene_camera
c
      subroutine compile_tria_chords(maxchords,maxcompile,
     & ntrii,nrknot,nx,ny,
     & xtrian,ytrian,necke,ixtri,iytri,crx,cry,bb,
     1 inpts,indchd,indit,indscl,filename,dataname,
     2 nxcam,xmin,xmax,nycam,ymin,ymax)

      use b2mod_types
      implicit none

      integer maxchords,maxcompile,ntrii,nrknot,nx,ny,inpts
      integer necke(3,ntrii),ixtri(ntrii),iytri(ntrii)
      integer indchd(maxchords+1),indit(maxcompile)
      real (kind=R8) ::
     * indscl(maxcompile),
     * xtrian(nrknot),
     * ytrian(nrknot),
     1 crx(-1:nx,-1:ny,0:3),
     2 cry(-1:nx,-1:ny,0:3),
     3 bb(-1:nx,-1:ny,0:3)
      character*(*) filename,dataname
      integer nxcam,nycam
      real (kind=R8) ::
     *   xmin,ymin,xmax,ymax

      logical open_file,vector
      integer ind,indold,npoints,lnblnk
      external lnblnk, subini, subend, xertst
      real (kind=R8) ::
     *   x1,x2,y1,y2,z1,z2

      call subini ('compile_tria_chords')
      if(.not.open_file(99,filename,
     1 'data.local/chords|data/chords')) then
        call subend ()
        return
      end if
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
      call chord_tria_compile(npoints,ntrii,nrknot,nx,ny,
     * xtrian,ytrian,necke,ixtri,iytri,crx,cry,bb,
     1 x1,y1,z1,x2,y2,z2,vector,
     1 indchd(inpts),indit(ind),indscl(ind))
      indold=ind
      ind=ind+indchd(inpts)
      call xertst(ind .le. maxcompile,
     1  'increase maxcompile in process_eirene_camera')
      indchd(inpts)=indold
      if(mod(inpts,1000) .eq. 0) 
     > write(*,*) 'Processed ',inpts,' chords'
      if(inpts .lt. maxchords) goto 10
      call xerrab ('increase array size in process_eirene_chords')
   30 write(*,*) 'Error in reading ',filename(1:lnblnk(filename))
      call subend ()
      return
   40 continue
      indchd(inpts+1)=ind
      write(*,*) ind-1,' elements of ',maxcompile,' used'
      close(99)
      if(inpts .ne. nxcam*nycam) then
        write(*,*) 'Problem with chords, ',nxcam*nycam,
     1   ' expected, ',inpts,' read in'
        call subend ()
        return
      endif
      call subend ()
      return
      end subroutine compile_tria_chords
c
#ifdef B25_EIRENE
      subroutine calculate_tzvalues(maxchords,maxcompile,
     1 inpts,indchd,indit,indscl,workarr,
     2 zvalues,filename,dataname,
     3 nxcam,xmin,xmax,nycam,ymin,ymax)

      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot
      use eirmod_ctrig
      implicit none

      integer maxchords,maxcompile,inpts
      integer indchd(maxchords+1),indit(maxcompile)
      real (kind=R8) ::
     * indscl(maxcompile)
      real (kind=R8) ::
     *   zvalues(maxchords)
      real (kind=R8) ::
     *   workarr(ntrii)
      character*(*) filename,dataname
      integer nxcam,nycam
      real (kind=R8) ::
     *   xmin,ymin,xmax,ymax

      logical open_file,vector
      integer npoints,lnblnk,i,j
      real (kind=R8) ::
     *   x1,x2,y1,y2,z1,z2,chord_tria_integral
      external lnblnk,chord_tria_integral, subini, subend, xerrab

      call subini ('calculate_tzvalues')
      if(use_compile) then
        do i=1,inpts
          zvalues(i)=0.0_R8
          do j=indchd(i),indchd(i+1)-1
            zvalues(i)=zvalues(i)+workarr(indit(j))*indscl(j)
          enddo
        enddo
      else
        if(.not.open_file(99,filename,
     1   'data.local/chords|data/chords')) then
          call subend ()
          return
        end if
        read(99,*,err=30,end=30) dataname,npoints,
     1   nxcam,xmin,xmax,nycam,ymin,ymax
        vector=npoints .lt. 0
        if(vector)
     1   dataname=dataname(1:lnblnk(dataname)) // ' (vector)'
        npoints=abs(npoints)
        inpts=0
   10   read(99,*,err=30,end=40) x1,y1,z1,x2,y2,z2
        inpts=inpts+1
        zvalues(inpts)=chord_tria_integral(ntrii,nrknot,nx,ny,workarr,
     1   xtrian,ytrian,necke,ixtri,iytri,crx,cry,bb,npoints,
     2   x1,y1,z1,x2,y2,z2,vector)
        if(mod(inpts,1000) .eq. 0)
     >   write(*,*) 'Processed ',inpts,' chords'
        if(inpts .lt. maxchords) goto 10
        call xerrab ('increase array size in process_eirene_chords')
   30   write(*,*) 'Error in reading ',filename(1:lnblnk(filename))
        call subend ()
        return
   40   continue
        close(99)
        if(inpts .ne. nxcam*nycam) then
          write(*,*) 'Problem with chords, ',nxcam*nycam,
     1     ' expected, ',inpts,' read in'
          call subend ()
          return
        endif
      endif
      call subend ()
      return
      end subroutine calculate_tzvalues
#endif
