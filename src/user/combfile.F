      subroutine combfile(lun,fn,lnfn,chdr,mhdr,nhdr,label,keep,xtrn,
     ,      cdfl,vdfl,ndfl,lwfs,lwfr,csep,mrks,nmrk,vtrl,ntrl,frmd,iret)
c
c  version : 20.02.2010 20:04
c
c=======================================================================
c*** Checks whether the specified tracing file exists and conforms to
c*** the current file format.
c*** If OK, then leaves the file open on the unit lun, positioned for
c*** append. If it does not exist, then opens a new file and writes the
c*** header there. If it does not conform, then tries to re-arrange the
c*** columns of the existing file, re-writes the file under a new name,
c*** padding the mising columns, and leaves it open on the unit lun,
c*** ready for append.
c*** This subroutine returns the actual file name, the number and values
c*** of the markers found in the source file, and the number and
c*** values of the trailing values to be used for padding the columns
c*** which are no longer present in the file.
c***
c*** The general format of a tracing file looks like the following:
c***  (a) A header consisting of 1 or more lines starting with "#". The
c***      last line of the header contains the list of column headers
c***      separated by white space. If the header consists of more than
c***      1 line, then the first line contains a list of markers -
c***      that is, integers identifying the data (e.g., version number,
c***      number of species, and so on). The identifier field may be
c***      limited by a second "#".
c***      All the intermediate header lines are ignored.
c***  (b) The data columns separated by white space. All the data are
c***      supposed to be in "e" format, the width of the first column
c***      may differ from the rest.
c***
c*** Input:
c***  lun   unit number to connect the file
c***  fn    file name (may be modified on output)
c***  lnfn  its length as declared in the calling routine
c***  chdr  array of column headers in actual format
c***  mhdr  maximum width of a column header
c***  nhdr  actual number of columns
c***  label a text label to appear in the output file
c***  keep  if true, then keep the obsolete columns as the trailing ones
c***  xtrn  external subroutine to be called to transform some old data
c***            call xtrn(data,chdr,mhdr,nhdr,mrks,nmrk,first)
c***          where data and chdr are for the re-arranged data row and
c***          column header, mrks are the markers found in the file,
c***          first=.true. for the first call.
c***          A dummy with no action is called "xtrn_dum"
c***  mrks  array of markers to be written in the file
c***  nmrk  actual number of markers
c***  cdfl  column names for nonzero defaults
c***  vdfl  corresponding default values
c***  ndfl  number of these values
c***        - default values pad the newly introduced columns in the
c***          old part of the file
c***  ntrl  maximum number of trailing columns
c***  lwfs  number of decimal digits in the first column (d in eW.d)
c***  lwfr  number of decimal digits in other columns
c***  csep  the column separator (space or tab)
c***
c*** Output:
c***  fn    actual file name
c***  mrks  array of markers found in the file
c***  nmrk  number of them in the file
c***  vtrl  array of values in the trailing columns
c***  ntrl  their actual number
c***        - trailing values pad the columns from the old file which
c***          are not present in the actual format
c***  frmd  format to be used for the data output (separators included)
c***  iret  return code:
c***           0 - the same file format or no file present (old action)
c***          <0 - unable to prepare the file
c***          >0 - the file changed. The rightmost bits (mask) mean:
c***                 1  columns added
c***                 2  columns removed (to be padded)
c***                 4  version (the first marker) changed
c***                 8  other changes in the markers
c***                16  new file opened
c=======================================================================
      use b2mod_types
      implicit none
      integer lun,lnfn,nhdr,mhdr,ndfl,nmrk,ntrl,iret,lwfs,lwfr,
     ,  mrks(nmrk)
      logical keep
      real (kind=R8) :: vdfl(ndfl),vtrl(ntrl)
      character*(mhdr) chdr(nhdr),cdfl(ndfl),fn*(*),label*(*),csep*1,
     ,  frmd*(*)
      external xtrn

c*** lsuf: number of digits in the numerical part of the file suffix
c*** nummax: maximum width of a numerical field in the file

      integer lsuf,nummax
      parameter (lsuf=2,nummax=24)
      logical ex,op,first
      character fbn*(lnfn),frmi22*8,chdf(2*nhdr)*(mhdr)
      character cbln*(nummax) / ' ' /
      character c,tab,f*(lnfn+3+lsuf),cc*20,frmhdr*12,frmdat*48
      integer i,j,k,l,m,n,ir,iv,kb,kt,lrd,nhd,nmr,ktrl,ltrl,
     ,        lnrw,lvrw,ncln,nrec,ivmax,mrkr(nmrk)
      integer lfrmd ! minimum length of the frmd variable
      parameter(lfrmd=36)
      integer ins_hlp(nhdr),ind_hlp(nhdr)
      integer,allocatable :: inr_hlp(:),int_hlp(:)
      character*(mhdr),allocatable :: chdw(:)
      character*(nhdr*48),allocatable :: cline(:)
      real (kind=R8), allocatable :: rowr(:),roww(:)
      logical streql
      external streql
c=======================================================================
c*** Prepare necessary constants

c!###{
c      write(0,*) 'combfile: lun,lnfn,mhdr,nhdr,keep,ndfl,lwfs,lwfr,',
c     ,  'nmrk,ntrl=',lun,lnfn,mhdr,nhdr,keep,ndfl,lwfs,lwfr,nmrk,ntrl
c!###}
      if(len(frmd).lt.lfrmd) then !{
        print *,'combfile: insufficient length of the frmd character ',
     ,                                    'variable - ',len(frmd),lfrmd
        iret = -4
        return
      end if !}
      tab=char(9)
      if(lsuf.le.0 .or. lsuf.gt.9) then !{
        print *,'combfile: invalid value of lsuf=',lsuf
        iret = -2
        return
      end if !}
      write(frmi22,'(4h(a,i,i1,1h.,i1,1h))') lsuf,lsuf
      i=1
      ivmax=0
      do j=1,lsuf !{
        ivmax=ivmax+i*9
        i=10*i
      end do !}
      lnrw=nhdr

c*** Check whether the file exists and find its actual name.
c*** When the routine has to re-arrange the file, it gives the file a
c*** new name, adding a suffix ".n.NN", where NN is a two-digit integer,
c*** to the base name (that is, without such suffix). The actual name is
c*** therefore either the composed name with the highest consecutive NN
c*** starting from 01, or the base name itself if no file "<fn>.n.01"
c*** is present.

      fbn=fn
      iv=0
      j=index(fn,'.',back=.true.)
      if(j.gt.3 .and. fn(j-2:j).eq.'.n.' .and.
     .                                  len_trim(fn).eq.j+lsuf) then !{
        read(fn(j+1:),'(i4)',err=10) i
        fbn=fn(1:j-3)
        iv=i
      end if !}
 10   continue

      ex=.true.
      inquire(file=fn,exist=ex)
c      write(0,*) 'ex,fn,iv,ivmax=',ex,trim(fn),iv,ivmax !###
      if(ex) then !{
        op=.true.
        i=iv
        do while (op .and. i.le.ivmax) !{
          i=i+1
          write(f,frmi22) trim(fbn)//'.n.',i
          inquire(file=f,exist=op)
c          write(0,*) 'op,f=',op,'  ',trim(f) !###
          if(op) fn=f
        end do !}
        if(i.gt.ivmax) then !{
          write(*,*) '!!! combfile: too many versions of ',trim(fbn)
          iret=-1
          return
        end if !}
        iv=i-1
c        write(*,*) 'combfile: reading ',trim(fn),' =>',iv !###
      end if !}

      if(ex) then !{

c*** Parse the header

        ir=-999
        do while (ir.ne.0 .and. iv.ge.0) !{
          ir=-999
c          write(0,*) 'to open: ',lun,': ',trim(fn) !###
          open(lun,file=fn,err=20)
c          write(0,*) 'opened: ',trim(fn) !###
          m=mhdr
          nhd=2*nhdr
          nmr=nmrk
          mrkr=mrks
          call parsehdr(lun,chdf,nhd,m,mrkr,nmr,ir)
c          write(0,*) 'after parsehdr. nhd,m,nmr,ir= ',nhd,m,nmr,ir !###
c          write(0,*) 'mrkr: ',mrkr !###
c          write(0,*) 'mrks: ',mrks !###
 20       if(ir.ne.0) then !{
            write(*,*) '!!!! ',trim(fn),' :',ir,' :',nhd,2*nhdr,' :',
     ,                                              m,mhdr,' :',nmr,nmrk
            if(ir.gt.0) then !{
              write(*,*) 'combfile: limits violation in the header of ',
     ,          trim(fn)
              do while (ex .and. iv.lt.ivmax) !{
                iv=iv+1
                write(fn,frmi22) trim(fbn)//'.n.',iv
                inquire(file=fn,exist=ex)
              end do !}
              if(ex) then !{
                write(*,*) 'combfile: too many versions of ',trim(fbn),
     ,            '  -  giving up!'
                iret=-1
                return
              end if !}
            else !}{
              if(iv.le.0) then !{
                iv=0
                ex=.false.
                ir=0
              else !}{
                iv=iv-1
              end if !}
c              write(0,*) 'iv=',iv  !###
              if(iv.gt.0) then !{
                write(fn,frmi22) trim(fbn)//'.n.',iv
              else !}{
                fn=fbn
              end if !}
            end if !}
            if(.not.ex)
     -        write(*,*) 'combfile: cannot continue - starting ',trim(f)
          end if !}
c          write(0,*) 'combfile: ir-iv loop. ir,iv=',ir,iv !###
        end do !}
      end if !}
c      write(0,*) 'ex,fn=',ex,'  ',trim(fn)  !###

      if(ex) then !{

c*** Check  whether the number of markers or their values have changed

        iret=0
        op=nmr.ne.nmrk
        if(op .and. (nmr.eq.0 .or. nmrk.eq.0) .or. nmr.gt.0 .and.
     .                                      mrks(1).ne.mrkr(1)) then !{
          iret=iret+4
          write(*,*) '!!! combfile: current version ',mrks(1),
     ,                         '  differs from that in the file ',fn
        end if !}
        if(.not.op) then !{
          do i=2,nmr !{
            if(mrks(i).ne.mrkr(i) .and.
     .                      (mrks(i).ne.0 .or. mrkr(i).ne.0)) op=.true.
          end do !}
        end if !}
        if(op) iret=iret+8
        if(op) write(*,*) '!!! combfile: markers change',nmrk,nmr

c*** Compare the headers
c*** At this point, the array chdr contains nhdr actual column headers
c*** and chdf contains nhd column headers from the file.

c*** ins_hlp for the indices of new columns in the old header,
c*** inr_hlp for the indices of old columns in the new header,
c*** int_hlp for the indices of trailing columns in the old header,
c*** ind_hlp for the indices of new columns in the padding set

        allocate(inr_hlp(nhd))
        allocate(int_hlp(nhd))

        ins_hlp=0
        inr_hlp=0
        int_hlp=0
        ind_hlp=0

        do i=1,nhdr !{
          do j=1,nhd !{
            if(chdr(i).eq.chdf(j)) then !{
              if(ins_hlp(i).eq.0) ins_hlp(i)=j
              if(inr_hlp(j).eq.0) inr_hlp(j)=i
            end if !}
          end do !}
        end do !}

c*** Defaults for padding

        op=.false.
        do i=1,nhdr !{
          if(ins_hlp(i).eq.0) then !{
            do j=1,ndfl !{
              if(chdr(i).eq.cdfl(j) .and. ind_hlp(i).eq.0) then !{
                ind_hlp(i)=j
                op=.true.
              end if !}
            end do !}
          end if !}
        end do !}
        if(op) iret=iret+1

c*** Trailing columns

        ktrl=0
        do j=1,nhd !{
          if(inr_hlp(j).eq.0) then !{
            ktrl=ktrl+1
            int_hlp(ktrl)=j
          end if !}
        end do !}
        if(ktrl.gt.0) iret=iret+2
        if(keep .and. ktrl.gt.ntrl) then !{
          print *,'combfile: too many trailing columns ',ktrl,ntrl

c***    look for the name for a new file

          do while (ex .and. iv.lt.ivmax) !{
            iv=iv+1
            write(fn,frmi22) trim(fbn)//'.n.',iv
            inquire(file=fn,exist=ex)
          end do !}
          if(ex) then !{
            write(*,*) 'combfile: too many versions of ',trim(fbn),
     ,        '  -  giving up!'
            iret=-2
            deallocate(inr_hlp,int_hlp)
            return
          end if !}
          iret=iret+16
        end if !}
        ltrl=ktrl
        if(.not.keep) ktrl=0
      end if !}
cc!###{
c      print *,'ex,ktrl=',ex,ktrl
c      print '(a/(20i4))','ins_hlp',ins_hlp
c      print '(a/(20i4))','inr_hlp',inr_hlp
c      print '(a/(20i4))','ind_hlp',ind_hlp
c      if(ktrl.gt.0) print '(a/(20i4))','int_hlp',(int_hlp(i),i=1,ktrl)
c      write(0,*) 'ex,iret=',ex,iret
cc!###}

      if(ex) then !{
        if(iret.gt.0) then !{

c*** Find a free i/o unit for re-writing the file

          op=.true.
          lrd=0
          do while (op .and. lrd.lt.98) !{
            lrd=lrd+1
            inquire(lrd,opened=op)
          end do !}
          do while (ex .and. iv.lt.ivmax) !{
            iv=iv+1
            write(f,frmi22) trim(fbn)//'.n.',iv
            inquire(file=f,exist=ex)
          end do !}
          if(ex) then !{
            write(*,*) 'combfile: too many versions of ',trim(fbn),
     ,        '  -  giving up!'
            iret=-3
            deallocate(inr_hlp,int_hlp)
            return
          end if !}
          close(lun)
          fn=f
          write(0,*) 'open lun',lun,trim(fn)
          open(lun,file=fn)
          write(0,*) 'open lrd',lrd,trim(fbn)
          open(lrd,file=fbn)
          rewind(lrd)
          c='#'
c          i=0 !###
          do while (c.eq.'#') !{
            read(lrd,'(a)',end=30) cc
            cc=adjustl(cc)
            c=cc(1:1)
c            i=i+1 !###
          end do !}
 30       backspace lrd
c          write(0,*) 'read comment lines: ',i  !###

          lnrw=nhdr+ktrl
          lvrw=nhdr+ltrl
          allocate(chdw(lvrw))
          allocate(roww(lvrw))
          allocate(rowr(nhd))

c*** Produce the header

          do i=1,nhdr !{
            chdw(i)=chdr(i)
          end do !}
          do i=1,ltrl !{
            chdw(i+nhdr)=chdf(int_hlp(i))
          end do !}
c          write(0,*) 'chdw',lnrw,' :',('  ',trim(chdw(i)),i=1,lvrw)!###
          if(nmrk.gt.0) then !{
            write(frmhdr,'(6h(2h# ,,i2,4hi12))') nmrk
            write(lun,frmhdr) mrks
            write(lun,'(2a)') '##  ',label
          end if !}
          write(frmhdr,'(5h(1h#,,i4,2ha))') 2*lnrw
          write(lun,frmhdr) (trim(chdw(i)),tab,i=1,lnrw)

c*** re-write the existing file

          ncln=1
          allocate(cline(ncln))
          write(frmhdr,'(1h(,i2,2ha))') ncln
          l=len(cline)
          nrec=0
          first=.true.
          write(frmdat,700) lwfs+8,lwfs,csep,lnrw-1,lwfr+8,lwfr,csep
c          print *,trim(frmdat)  !###
          do while (.true.) !{
            read(lrd,frmhdr,end=45) cline
            if(.not.streql(cline(ncln)(l-nummax:l),' ')) then !{
              if(ncln.gt.1) then !{

c*** Check whether any data fields are split between the line parts
c*** and fix this if possible

                do j=ncln,2,-1 !{
                  if(cline(j)(1:1).ne.' ' .and.
     .                      cline(j)(1:1).ne.tab .and.
     .                            cline(j-1)(l:l).ne.' ' .and.
     .                                  cline(j-1)(l:l).ne.tab) then !{
c*** non-space on both sides of the section borders - the field is split
                    kb=index(cline(j-1),' ',back=.true.)
                    kt=index(cline(j-1),tab,back=.true.)
                    n=l-max(kb,kt)  ! length of the 1st field fragment
                    k=index(cline(j),cbln)
                    kt=k+nummax
                    kb=k+1
                    if(k.gt.0) then !{
                      if(.not.streql(cline(j)(k:l),' ')) then !{ non-empty tail
                        cline(j)(k:k)=csep
                        do while ( k.gt.0 .and.
     .                            .not.streql(cline(j)(k:l),' ')) !{
                          cline(j)(kb:kb+l-kt)=cline(j)(kt:l)
                          k=index(cline(j),cbln)
                          kb=k
                          kt=k+nummax
                        end do !}
                      end if !}
                      if(.not.streql(cline(j)(k:l),' ')) then !{ shift to the right
                        cline(j)(n+1:l)=cline(j)(1:l-n)
                        cline(j)(1:n)=cline(j-1)(l-n+1:l)
                        cline(j-1)(l-n+1:l)=repeat(' ',n)
                      end if !}
                    end if !}
                    if(k.eq.0) then !{ not enough space -> next section
                      print *,'combfile: next section not ready'  !###
                      print *,'nrec,l,ncln=',nrec,l,ncln          !###
                      stop                                        !###
                    end if !}
                  end if !}
                end do !}
              end if !}
              cline(ncln)=trim(cline(ncln))//'/'
              read(cline,*,err=40) rowr
              nrec=nrec+1
              do i=1,nhdr !{
                j=ins_hlp(i)
                if(j.gt.0) then !{
                  roww(i)=rowr(j)
                else if(ind_hlp(i).gt.0) then !}{
                  roww(i)=vdfl(ind_hlp(i))
                else !}{
                  roww(i)=0.0_R8
                end if !}
              end do !}
              do i=1,ltrl !{
                roww(i+nhdr)=rowr(int_hlp(i))
              end do !}
              call xtrn(roww,chdw,mhdr,lvrw,mrkr,nmr,first)
              first=.false.
              write(lun,frmdat) (roww(i),i=1,lnrw)
              nrec=nrec+1
            else !}{ Probably, the line is too short - add a section
              ncln=ncln+1
              deallocate(cline)
              allocate(cline(ncln))
              write(frmhdr,'(1h(,i2,2ha))') ncln
              backspace lrd
            end if !}
 40       continue
          end do !}
 45       continue
          close(lrd)
          if(lnrw.eq.lvrw) then !{
            ntrl=lnrw-nhdr
            do i=1,ntrl !{
              vtrl(i)=roww(i+nhdr)
            end do !}
          else !}{
            ntrl=0
          end if !}
        else !}{

c*** Read to the end of the file

 50      read(lun,'(a1)',end=55) c
         go to 50
 55      backspace lun
        end if !}

      else !}{

c*** Open a new tracing file

        write(0,*) 'opening the tracing file ',trim(fn)  !###
        open(lun,file=fn)

c*** Write the header

        if(nmrk.gt.0) then !{
          write(frmhdr,'(6h(2h# ,,i2,4hi12))') nmrk
          write(lun,frmhdr) mrks
          write(lun,'(2a)') '##  ',label
        end if !}
        write(frmhdr,'(5h(1h#,,i4,2ha))') 2*nhdr
        write(lun,frmhdr) (trim(chdr(i)),tab,i=1,nhdr)
        ntrl=0
        iret=0
      end if !}

      if(allocated(inr_hlp)) deallocate(inr_hlp,int_hlp)
      if(allocated(chdw)) deallocate(chdw,roww,rowr)
      if(allocated(cline)) deallocate(cline)

      write(frmd,700)lwfs+8,lwfs,csep,lnrw-1,lwfr+8,lwfr,csep
  700 format(8h(1p,1x,e,i2,1h.,i2,3h,1h,a1,1h,,
     ,  i4,2h(e,i2,1h.,i2,3h,1h,a1,2h)))
c      write(0,*) 'leaving combfile. lnrw,frmd=',lnrw,'  ',trim(frmd) !###
c=======================================================================
      end


c*** Dummy data transformation routine

      subroutine xtrn_dum(data,chdr,mhdr,nhdr,mrks,nmrk,first)
      use b2mod_types
      implicit none
      real (kind=R8) :: data(nhdr)
      integer :: nmrk, nhdr, mhdr
      integer :: mrks(nmrk)
      character*(mhdr) chdr(nhdr)
      logical first
      return
      end
