c
      subroutine process_chord_slice(filename)
      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot
      implicit none
      character*(*) filename
      character*80 dataname,tmpstring
      character*256 extralabel
      logical open_file,vector
      integer lngarr
      parameter (lngarr=500 000)
      real (kind=R8) ::
     *   xvalues(lngarr),yvalues(lngarr),chord_value,
     1   xmin,xmax,ymin,ymax,x1,y1,z1,dx1,dy1,dz1
      integer inpts,i,is,npoints,lnblnk
      external lnblnk,chord_value

      do is=0,ndim(3,iwrk)-1
        if(.not.open_file(99,filename,
     1   'data.local/chords|data/chords')) return
        read(99,*,err=30,end=30) dataname,npoints
        vector=npoints .lt. 0
        dx1=0
        dy1=0
        dz1=0
        inpts=0
   10   if(vector) then
          read(99,*,err=30,end=40) x1,y1,z1,dx1,dy1,dz1,xvalues(inpts+1)
        else
          read(99,*,err=30,end=40) x1,y1,z1,xvalues(inpts+1)
        endif
        inpts=inpts+1
        yvalues(inpts)=chord_value(nx,ny,
     1   work(-1,-1,is,iwrk),crx,cry,bb,
     2   x1,y1,z1,
     3   dx1,dy1,dz1,vector)
        if(inpts .lt. lngarr) goto 10
        stop 'increase array size in process_chord_slice'
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

        tmpstring=textwk(is,iwrk)

        call init_page(
     1   graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk))) //
     2   tmpstring(1:lnblnk(tmpstring)) // ' '//
     2   dataname(1:lnblnk(dataname)),
     1   xmin,xmax,ymin,ymax,linlog,.false.,.false.,.false.)
        write(*,*) 'xmin,xmax,ymin,ymax ',xmin,xmax,ymin,ymax
        if(writing_to_file) then
          write(nwtf,'(''# '',a)')                   !srv 21.09.11
     1     graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk))) //
     2     tmpstring(1:lnblnk(tmpstring)) // ' '//
     2     dataname(1:lnblnk(dataname))
          call get_extralabel(extralabel)
          if(len_trim(extralabel).ne.0) then
            write(nwtf,'(''# '',a)') trim(extralabel)                   !srv 21.09.11
          endif
          do i=1,inpts
            write(nwtf,*) xvalues(i),yvalues(i)                   !srv 21.09.11
          enddo
        endif
        do i=1,inpts-1
          call lined_d(xvalues(i),yvalues(i),xvalues(i+1),yvalues(i+1))
        enddo
        call end_page
        close(99)
      enddo
      return
      end subroutine process_chord_slice
c
      subroutine do_average(val,no)
      use b2mod_types
      use b2mod_b2plot
      implicit none
      integer no
      real (kind=R8) ::
     *   val(no)
      integer i,j
      real (kind=R8) ::
     *   ave

      do i=1,no-average+1,average
        ave=0
        do j=i,i+average-1
          ave=ave+val(j)
        enddo
        val(i/average+1)=ave/average
      enddo
      return
      end subroutine do_average
