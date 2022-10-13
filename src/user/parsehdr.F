c      integer nn,mm,nm
cc      parameter (nn=200,mm=16, nm=12)
c      parameter (nn=100,mm=10, nm=2)
c      integer mr(nm)
c      character chdr(nn)*(mm)
c      open(1,file='integral.trc')
cc      open(1,file='blnn.trc')
cc      open(1,file='traces')
c      n=nm
c      ll=nn
c      kk=mm
c      call parsehdr(1,chdr,ll,mr,n,ir)
c      print *,'lhdr,khdr,nm,iret=',ll,kk,n,ir
c      if(nm.gt.0) print *,(mr(i),i=1,n)
c      if(ir.eq.0 .and. ll.gt.0) print '(1x,a)',(trim(chdr(i)),i=1,ll)
c      end



      subroutine parsehdr(lun,chdr,nhdr,mrks,nmrk,iret)
c
c  version : 22.11.2006 00:41
c
c=======================================================================
c*** Parses the header lines of a file pre-connected to the unit lun
c*** The header lines starrt with a "#",
c*** the first word of the first line is the version number,
c*** the last line contains the column headers
c***
c*** Input:
c***
c***  lun   the unit number where the file is pre-connected
c***  nhdr  the dimension of the chdr array in the calling routine
c***  mhdr  the length of an element of chdr
c***  nmrk  the maximum number of file markers
c***
c*** Output:
c***
c***  chdr  the column headers
c***  nhdr  the number of column headers
c***  mrks  the array of file markers
c***  nmrk  the actual number of file markers
c***  iret  the return code: 0 - OK,
c***    -1 - no file, -2 - binary file, -3 - general file error,
c***    -4 - bad header, -5 - too long header line.
c***    In case of iret > 0, the rightmost bits (mask) mean:
c***    1 - too many columns, 2 - too wide headers, 4 - too many markers
c=======================================================================
      use b2mod_diag
      implicit none
      integer lun,nhdr,nmrk,iret,mrks(nmrk)
      character*(mhdr) chdr(nhdr)
      character c,tab,line*(nhdr*(mhdr+1)*2),cline*(nhdr*(mhdr+1)*2+10)
      integer i,k,l,m,n
      logical lopn
      character*3 cfrm,cseq
      logical streql
      intrinsic trim,len_trim,adjustl,huge
      external streql
c=======================================================================

      tab=char(9)
c      print *,'parsehdr: lun,nhdr,mhdr,len=',lun,nhdr,mhdr,len(line)!###

c*** Check the file
      iret=0
      inquire(lun,opened=lopn,sequential=cseq,formatted=cfrm,err=999)
      if(.not.lopn) iret=-1
      if(.not.streql(cseq,'YES')) iret=-2
      if(.not.streql(cfrm,'YES')) iret=-2
      if(iret.ne.0) return

c*** The file is opened and of right format. Check the header

      rewind(lun,err=999)
      read(lun,'(a)',end=998,err=999) line
      line=adjustl(line)
      if(line(1:1).ne.'#') then !{
        iret=-4
        return
      end if !}
      line(1:1)=' '
      line=adjustl(line)
      l=index(line,'#')
      if(l.gt.0) line=line(1:l-1)
      if(l.le.0) l=len_trim(line)
      if(l.le.0) go to 998
      cline=line

      iret=0
      m=-huge(l)
      mrks=-huge(l)
      k=0
      c='#'
      do while (c.eq.'#') !{
        read(lun,'(a)',end=10,err=999) line
        k=k+1
        line=adjustl(line)
        c=line(1:1)
      end do !}
   10 if(k.gt.1) then !{
        cline=trim(cline)//'/'
        read(cline,*,err=15) mrks,m
        if(m.ne.-huge(l)) iret=iret+4
        m=nmrk
c        do while (m.lt.nmrk .and. mrks(m+1).ne.-huge(l)) !{}
        do i=1,nmrk !{
          if(mrks(i).eq.-huge(l)) then !{
            m=i-1
            exit
          end if !}
        end do !}
        nmrk=m
      else if(k.le.0) then !}{
        go to 998
      else !}{
        nmrk=0
      end if !}
      go to 20
   15 nmrk=0
   20 backspace(lun,err=999)
      backspace(lun,err=999)
      read(lun,'(a)',err=998) cline
      line=cline
      if(.not.streql(line,cline)) then !{
        iret=-5
        return
      end if !}
      line=adjustl(line)
      if(line(1:1).ne.'#') then !{
        iret=-4
        return
      end if !}
      n=0
      k=0
      line(1:1)=' '
      do while (.not.streql(line,' ')) !{
        line=adjustl(line)
        if(line(1:1).eq.tab) then !{
          line(1:1)=' '
        else !}{
          l=index(line,tab)
          m=index(line,' ')
          if(l.gt.0 .and. m.gt.0) then !{
            l=min(index(line,tab),index(line,' '))-1
          else !}{
            l=max(index(line,tab),index(line,' '))-1
          end if !}
          if(l.le.0) l=len_trim(line)
          k=max(l,k)
          n=n+1
          if(n.le.nhdr) chdr(n)=line(1:l)
          line=line(l+2:len_trim(line))
        end if !}
      end do !}
      if(k.gt.mhdr) iret=iret+2
      if(n.gt.nhdr) iret=iret+1
      nhdr=n
!xpb      mhdr=k ! This was overwriting a parameter!
      return
c=======================================================================
  998 iret=-4
      return
  999 iret=-3
      return
c=======================================================================
      end subroutine parsehdr
