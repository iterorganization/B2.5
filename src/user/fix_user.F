      subroutine fix_user(data,chdr,nhdr,mrks,nmrk,first)

c  version : 04.04.2006 16:59

c=======================================================================
c*** A filter to be used with the combfile subroutine to fix the
c*** user.trc file when the b2usrtrc subroutine changes.
c*** Input
c***  data  the data row, already amended by combfile
c***  chdr  the data header, already amended by combfile
c***  nhdr  the number of positions in the data and header rows
c***  mrks  the markers as found in the file
c***  nmrk  the number of markers found in the file
c***  first if true, then this is the first call
c-----------------------------------------------------------------------

      use b2mod_types
      use b2mod_constants
      use b2mod_diag
      implicit none
      integer nhdr,nmrk,mrks(nmrk)
      logical first
      real (kind=R8) :: data(nhdr)
      character*(mhdr) chdr(nhdr)
      logical streql
      external streql

c-----------------------------------------------------------------------
c*** lv     version number from the file (the first marker)
c*** lvn    starting from this version, no corrections
c*** nntr   number of diagnostic surface groups for neutral data
c*** number of neutral species found
c*** larea  positions of the surface area of the groups to be fixed
c*** lntntr positions of temperature of species i at surface group j
c*** lnnntr positions of density of species i at surface group j
c*** lprntr positions of pressure of species i at surface group j
c*** lspntr positions of pumping speed for species i at surface group j
c*** cspntr identifiers of the neutral species
c*** csrntr identifiers of the surface groups

      integer lv,lvn, nntr,nnsp
      integer,allocatable :: larea(:),lnarea(:),lntntr(:,:),lnnntr(:,:),
     ,  lprntr(:,:),lspntr(:,:)
      character*8,allocatable :: cspntr(:),csrntr(:)
      data lv /0/
      parameter (lvn=20060402)
      integer i,j,k,l,m,ii
      logical ex
c      character*48 frh,frd
      character*(mhdr) cc

      save
c=======================================================================
c*** prior to 20060402: fix the neutral pressure and pumping speed at
c***                    diagnostic surfaces, if their area is written
c***                    as "area", not "Area"

      if(lv.ge.20060402) return

      if (first) then !{
c!###{
c        print *
c        print *,'fix_user: before initialisation. ',
c     ,    'mhdr,nhdr,nmrk=',mhdr,nhdr,nmrk
c        if(nmrk.gt.0) then !{
c          print *,'mrks:',mrks
c        end if !}
c        print *
c        do i=1,nhdr !{
c          print '(1x,i3,2x,a,1p,e11.3)',i,chdr(i),data(i)
c        end do !}
c        print *
c!###}
        lv=mrks(1)
        first=.false.
        if(lv.ge.20060402) return

c*** Determine the number of surface groups to fix

        nntr=0
        k=0
        do i=1,nhdr !{
          if(index(chdr(i),'area_').eq.1) then !{
            nntr=nntr+1
          else if(index(chdr(i),'Area_').eq.1) then !}{
            k=k+1
          end if !}
        end do !}
        if(nntr.gt.0 .and. k.gt.0 .and. k.ne.nntr) then !{
          print *,'!!! fix_user: different number of Area and area ',
     ,                                                'headers ',k,nntr
        end if !}

c*** Determine the number of neutral species to fix

        if(nntr.gt.0) then !{
          nnsp=0
          do i=1,nhdr !{
            if(index(chdr(i),'area_').eq.1) then !{
              cc=chdr(i)(6:)
              l=len_trim(cc)
c              print *,'found :',chdr(i),':',cc,':',l !###
              k=0
              do j=1,nhdr !{
                if(index(chdr(j),'n_').eq.1) then !{
c!###{
c                  print *,'comparing :',chdr(j),':',trim('_'//cc),':',
c     ,              index(chdr(j),trim('_'//cc),back=.true.),
c     ,              len_trim(chdr(j))-l
c!###}
                  if(index(chdr(j),trim('_'//cc),back=.true.).eq.
     .                                      len_trim(chdr(j))-l) then !{
                    k=k+1
                  end if !}
                end if !}
              end do !}
              nnsp=max(nnsp,k)
            end if !}
          end do !}
c          print *,'fix_user (a): nntr,nnsp=',nntr,nnsp  !###
          if(nnsp.gt.0) then !{
            allocate(larea(nntr))
            allocate(lnarea(nntr))
            allocate(lntntr(nnsp,nntr))
            allocate(lnnntr(nnsp,nntr))
            allocate(lprntr(nnsp,nntr))
            allocate(lspntr(nnsp,nntr))
            allocate(cspntr(nnsp))
            allocate(csrntr(nntr))
            larea=0
            lnarea=0
            lntntr=0
            lnnntr=0
            lprntr=0
            lspntr=0
            cspntr=' '
            csrntr=' '
          else !}{
            print *,'!!! fix_user: no neutral species, but nntr=',nntr
            nntr=0
          end if !}
        end if !}
c        print *,'fix_user (b): nntr,nnsp=',nntr,nnsp  !###
        if(nntr.gt.0) then !{

c*** Fill the species names and find the positions of related quantities

          k=0
          do i=1,nhdr !{
            if(index(chdr(i),'area_').eq.1) then !{
              k=k+1
              csrntr(k)=chdr(i)(6:)
              larea(k)=i
            end if !}
          end do !}
          do i=1,nhdr !{
            if(index(chdr(i),'Area_').eq.1) then !{
              cc=chdr(i)(6:)
              do j=1,nntr !{
                if(streql(cc,csrntr(j))) then !{
                  lnarea(j)=i
                  exit
                end if !}
              end do !}
            end if !}
          end do !}
          k=0
          do i=1,nhdr !{
            if(index(chdr(i),'n_').eq.1 .or.
     .                index(chdr(i),'T_').eq.1 .or.
     .                      index(chdr(i),'p_').eq.1 .or.
     .                            index(chdr(i),'S_').eq.1) then !{
              m=index(chdr(i),'_')
              l=index(chdr(i)(m+1:),'_')+m
              cc=chdr(i)(l+1:)
              do ii=1,nntr !{
                ex=streql(cc,csrntr(ii))
                if(ex) exit
              end do !}
              if(ex) then !{
                cc=chdr(i)(m+1:l-1)
                do j=1,k !{
                  if(streql(cc,cspntr(j))) then !{
                    ex=.false.
                    exit
                  end if !}
                end do !}
                if(ex) then !{
                  k=k+1
                  j=k
                  cspntr(k)=cc
                end if !}
                if(chdr(i)(1:1).eq.'n') then !{
                  lnnntr(j,ii)=i
                else if(chdr(i)(1:1).eq.'T') then !}{
                  lntntr(j,ii)=i
                else if(chdr(i)(1:1).eq.'p') then !}{
                  lprntr(j,ii)=i
                else if(chdr(i)(1:1).eq.'S') then !}{
                  lspntr(j,ii)=i
                end if !}
              end if !}
            end if !}
          end do !}
        end if !}
c!###{
c        if(nntr.gt.0) then !{
c          write(frh,'(7h(1x,a8,,i2,3ha9))') nnsp
c          write(frd,'(10h(1x,a8,1p,,i2,3hi9))') nnsp
c          print frh,'lnnntr:',cspntr
c          do j=1,nntr !{
c            print frd,csrntr(j),(lnnntr(i,j),i=1,nnsp)
c          end do !}
c          print frh,'lntntr:',cspntr
c          do j=1,nntr !{
c            print frd,csrntr(j),(lntntr(i,j),i=1,nnsp)
c          end do !}
c          print frh,'lprntr:',cspntr
c          do j=1,nntr !{
c            print frd,csrntr(j),(lprntr(i,j),i=1,nnsp)
c          end do !}
c          print frh,'lspntr:',cspntr
c          do j=1,nntr !{
c            print frd,csrntr(j),(lspntr(i,j),i=1,nnsp)
c          end do !}
c        end if !}
c!###}
      end if !}
      if(nntr.gt.0) then !{

c*** Correct the neutral data

        do j=1,nntr !{
          m=larea(j)
          k=lnarea(j)
          if(m.gt.0 .and. k.gt.0) data(k)=data(m)
          do i=1,nnsp !{
            m=lprntr(i,j)
            k=lnnntr(i,j)
            l=lntntr(i,j)
            ii=lspntr(i,j)
            if(k.gt.0 .and. l.gt.0 .and. m.gt.0) then !{
              data(m)=ev*data(k)*data(l)
            end if !}
            if(ii.gt.0) data(ii)=1.283*data(ii)
          end do !}
        end do !}
      end if !}
c=======================================================================
      end
