c
      subroutine process_chords_boxcar(filename)
      use b2mod_types
      use b2us_data
      use b2mod_b2plot
      use b2mod_b2plot_debug
      use b2mod_b2plot_values
      use b2mod_subsys
      implicit none
      character*(*) filename
      character*80 dataname
      character*256 extralabel,tmpstring,titlestring
      logical open_file,vector
      integer lngarr
      parameter (lngarr=100 000)
      real (kind=R8) ::
     *   xvalues(lngarr),yvalue(lngarr),yvalues(lngarr),
     1   xmin,xmax,ymin,ymax,x1,y1,z1,x2,y2,z2,chord_integral
      integer inpts,i,is,npoints,lnblnk
      external lnblnk,chord_integral,xerrab

      call subini('process_chords_boxcar')
      do is=0,ndim(2,iwrk)-1
        if(.not.open_file(99,filename,
     1   'data.local/chords|data/chords')) then
          call subend ()
          return
        endif
        read(99,*,err=30,end=30) dataname,npoints
        vector=npoints .lt. 0
        if(vector) dataname=dataname(1:lnblnk(dataname)) // ' (vector)'
        npoints=abs(npoints)
        inpts=0
   10   read(99,*,err=30,end=40) x1,y1,z1,x2,y2,z2,xvalues(inpts+1)
        inpts=inpts+1
        yvalue(inpts)=chord_integral(mpg%nCv,work(1,is,iwrk),
     1   npoints,x1,y1,z1,x2,y2,z2,vector)
        if(inpts .lt. lngarr) goto 10
        call xerrab ('increase array size in process_chord')
   30   write(*,*) 'Error in reading ',filename(1:lnblnk(filename))
        call subend ()
        return
   40   continue
        close(99)
        if (inpts .lt. 2) write (*,*)
     .    ' Not enough chords for proper averaging. Data unreliable.'
        yvalues(1)=0.5_R8*(yvalue(1)+yvalue(2))
        yvalues(inpts)=0.5_R8*(yvalue(inpts)+yvalue(inpts-1))
        do i=2,inpts-1
          yvalues(i)=(yvalue(i-1)+yvalue(i)+yvalue(i+1))/3.0_R8
        enddo
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
        tmpstring=graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk))) //
     2   textwk(is,iwrk)(1:lnblnk(textwk(is,iwrk))) // ' ' //
     3   dataname(1:lnblnk(dataname)) // ' (Boxcar)'
        call customise_label(tmpstring,titlestring)
        call init_page(titlestring,
     4   xmin,xmax,ymin,ymax,linlog,.false.,.false.,.false.)
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
      return
      end subroutine process_chords_boxcar
