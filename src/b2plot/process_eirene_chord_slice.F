c
      subroutine process_eirene_chord_slice
#ifdef B25_EIRENE
     . (filename)
      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot
      use b2mod_b2plot_debug
      use b2mod_b2plot_values
      use b2mod_subsys
      use eirmod_ctrig
      use b2us_data
      implicit none
      character*(*) filename
      character*80 dataname
      character*256 extralabel,titlestring,tmpstring
      logical open_file,vector
      integer lngarr
      parameter (lngarr=500 000)
      real (kind=R8) ::
     *   xvalues(lngarr),yvalues(lngarr),chord_tria_value,
     1   xmin,xmax,ymin,ymax,x1,y1,z1,dx1,dy1,dz1
      integer inpts,i,is,npoints,lnblnk
      external lnblnk,chord_tria_value,xerrab

      call subini ('process_chord_slice')
      do is=0,nedim(2,iewrk)-1
        if(.not.open_file(99,filename,
     1   'data.local/chords|data/chords')) then
          call subend ()
          return
        end if
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
        yvalues(inpts)=chord_tria_value(geo,ntrii,nrknot,
     1   eirene_stack(1,is,iewrk),
     2   xtrian,ytrian,necke,ixtri,
     2   x1,y1,z1,dx1,dy1,dz1,vector)
        if(inpts .lt. lngarr) goto 10
        call xerrab
     .   ('increase array size in process_eirene_chord_slice')
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

        tmpstring=
     1   eirene_label(iewrk)(1:lnblnk(eirene_label(iewrk))) //
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
      end subroutine process_eirene_chord_slice
