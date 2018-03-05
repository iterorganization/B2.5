      subroutine xxdata_11 ( iunit  , iclass ,
     &                       isdimd , iddimd , itdimd ,
     &                       ndptnl , ndptn  , ndptnc , ndcnct ,
     &                       iz0    , is1min , is1max ,
     &                       nptnl  , nptn   , nptnc  ,
     &                       iptnla , iptna  , iptnca ,
     &                       ncnct  , icnctv ,
     &                       iblmx  , ismax  , dnr_ele, dnr_ams,
     &                       isppr  , ispbr  , isstgr ,
     &                       idmax  , itmax  ,
     &                       ddens  , dtev   , drcof  ,
     &                       lres   , lstan  , lptn
     &                     )
      use b2mod_types
      implicit none
c-----------------------------------------------------------------------
c
c  ****************** fortran77 subroutine: xxdata_11 ******************
c
c  purpose: to read a complete  adf11 file, check its class and
c           determine its standard, resolved and partition organisation.
c
c  calling program: various
c
c  notes:    (1) A 'standard' adf11 file contains gcr data between one
c                whole ionisation stage and another whole ionisation
c                stage.
c                A 'resolved' (or partial) adf11 file contains gcr data
c                between a set of metastables of one ionisation stage
c                and a set of metastables of another ionisation stage.
c                A resolved file is distinguished from a standard file
c                by the presence of a 'connection vector' in the adf11
c                data file header lines.
c                The connection vector specifies the number of meta-
c                stables in each ionisation stage which are coupled
c                together by gcr data.
c            (2) A 'partitioned' adf11 file contains gcr data between
c                clumps of ionisation stages or metastables or comb-
c                inations of the two called 'partitions'.
c                A 'partition level' is a specification of the
c                partitions which span all the ionisation stages (and
c                metastables) of an element.  Successive partition
c                levels give a hierarchy corresponding to larger
c                partitions and greater clumping.
c                A 'superstage' is a set of partitions which are close-
c                coupled.
c                There are thus equivalences :
c                        ionisation stage   -   superstage
c                        metastable         -   partition
c                        ion charge         -   superstage index
c                A partitioned adf11 file may be standard (with each
c                superstage comprising only one partition) or resolved.
c                A partitioned file is distinguished by the presence of
c                'partition specification block' in the adf11 data
c                file header lines.
c            (3) When a partition specification block is present, it
c                should be ordered from the highest partition level
c                index to lowest partition level index.  Thus the first
c                partition in the partition block has the least number
c                of partitions and the last has the greatest number.
c            (4) Twelve classes of adf11 data file may be read by the
c                subroutine as follow:
c
c                class index    type      GCR data content
c                -----------    ----      ----------------
c                    1          acd       recombination coeffts
c                    2          scd       ionisation coeffts
c                    3          ccd       CX recombination coeffts
c                    4          prb       recomb/brems power coeffts
c                    5          prc       CX power coeffts
c                    6          qcd       base meta. coupl. coeffts
c                    7          xcd       parent meta. coupl. coeffts
c                    8          plt       low level line power coeffts
c                    9          pls       represent. line power coefft
c                   10          zcd       effective charge
c                   11          ycd       effective squared charge
c                   12          ecd       effective ionisation potential
c
c            (5) A resolved adf11 file, with a connection vector, has a set
c                of names and pointers at precise positions in the data file
c                which are recognised.
c                The names are different for partitioned and unpartitioned
c                data files as follow:
c
c                       file      unpartitioned         partitioned
c                       class        names                 names
c
c                       (all)         z1                   s1
c
c                                (indices 1 and 2)     (indices 1 and 2)
c                       ----       ----     ----         ----     ----
c                       acd        iprt     igrd         ispp     ispb
c                       scd        iprt     igrd         ispp     ispb
c                       ccd        iprt     igrd         ispp     ispb
c                       prb        iprt                  ispp
c                       prc        iprt                  ispp
c                       qcd        igrd     jgrd         ispb     jspb
c                       xcd        iprt     jprt         ispp     jspp
c                       plt        igrd                  ispb
c                       pls        igrd                  ispb
c                       zcd        igrd                  ispb
c                       ycd        igrd                  ispb
c                       ecd        igrd                  ispb
c
c             (6) In partitioned nomenclature: s=superstage; p=partition;
c                 b=base (current superstage), p=parent (next up super-
c                 stage), c=child (next down superstage). Thus arrays
c                 `iprtr' and `igrd' in old notation are now substituted
c                 by `isppr' and `ispbr' respectively internally and in
c                 external naming.
c
c
c  subroutine:
c
c  input : (i*4)  iunit     = unit to which input file is allocated
c  input : (i*4)  iclass    = class of data (1 - 10 ):
c                               1-acd, 2-scd, 3-ccd, 4-prb, 5-prc
c                               6-qcd, 7-xcd, 8-plt, 9-pls,10-zcd
c                              11-ycd,12-ecd
c
c  input : (i*4)  isdimd    = maximum number of (sstage, parent, base)
c                             blocks in isonuclear master files
c  input : (i*4)  iddimd    = maximum number of dens values in
c                             isonuclear master files
c  input : (i*4)  itdimd    = maximum number of temp values in
c                             isonuclear master files
c  input : (i*4)  ndptnl    = maximum level of partitions
c  input : (i*4)  ndptn     = maximum no. of partitions in one level
c  input : (i*4)  ndptnc    = maximum no. of components in a partition
c  input : (i*4)  ndcnct    = maximum number of elements in connection
c                             vector
c
c  output: (i*4)  iz0       = nuclear charge
c  output: (i*4)  is1min    = minimum ion charge + 1
c                             (generalised to connection vector index)
c  output: (i*4)  is1max    = maximum ion charge + 1
c                             (note excludes the bare nucleus)
c                             (generalised to connection vector index
c                              and excludes last one which always remains
c                              the bare nucleus)
c  output: (i*4)  nptnl     = number of partition levels in block
c  output: (i*4)  nptn()    = number of partitions in partition level
c                             1st dim: partition level
c  output: (i*4)  nptnc(,)  = number of components in partition
c                             1st dim: partition level
c                             2nd dim: member partition in partition level
c  output: (i*4)  iptnla()  = partition level label (0=resolved root,1=
c                                                      unresolved root)
c                             1st dim: partition level index
c  output: (i*4)  iptna(,)  = partition member label (labelling starts at 0)
c                             1st dim: partition level index
c                             2nd dim: member partition index in partition
c                             level
c  output: (i*4)  iptnca(,,)= component label (labelling starts at 0)
c                             1st dim: partition level index
c                             2nd dim: member partition index in partition
c                             level
c                             3rd dim: component index of member partition
c  output: (i*4)  ncnct     = number of elements in connection vector
c  output: (i*4)  icnctv()  = connection vector of number of partitions
c                             of each superstage in resolved case
c                             including the bare nucleus
c                             1st dim: connection vector index
c
c  output: (i*4)  iblmx     = number of (sstage, parent, base)
c                             blocks in isonuclear master file
c  output: (i*4)  ismax     = number of charge states
c                             in isonuclear master file
c                             (generalises to number of elements in
c                              connection vector)
c  output: (c*12) dnr_ele   = CX donor element name for iclass = 3 or 5
c                             (blank if unset)
c  output: (r*8)  dnr_ams   = CX donor element mass for iclass = 3 or 5
c                             (0.0d0 if unset)
c  output: (i*4)  isppr()   = 1st (parent) index for each partition block
c                             1st dim: index of (sstage, parent, base)
c                                      block in isonuclear master file
c  output: (i*4)  ispbr()   = 2nd (base) index for each partition block
c                             1st dim: index of (sstage, parent, base)
c                                      block in isonuclear master file
c  output: (i*4)  isstgr()  = s1 for each resolved data block
c                             (generalises to connection vector index)
c                             1st dim: index of (sstage, parent, base)
c                                      block in isonuclear master file
c
c  output: (i*4)  idmax     = number of dens values in
c                             isonuclear master files
c  output: (i*4)  itmax     = number of temp values in
c                             isonuclear master files
c  output: (r*8)  ddens()   = log10(electron density(cm-3)) from adf11
c  output: (r*8)  dtev()    = log10(electron temperature (eV) from adf11
c  output: (r*8)  drcof(,,) = if(iclass <=9):
c                                log10(coll.-rad. coefft.) from
c                                isonuclear master file
c                             if(iclass >=10):
c                                coll.-rad. coefft. from
c                                isonuclear master file
c                             1st dim: index of (sstage, parent, base)
c                                      block in isonuclear master file
c                             2nd dim: electron temperature index
c                             3rd dim: electron density index

c  output: (l*4)  lres      = .true.  => partial file
c                           = .false. => not partial file
c  output: (l*4)  lstan     = .true.  => standard file
c                           = .false. => not standard file
c  output: (l*4)  lptn      = .true.  => partition block present
c                           = .false. => partition block not present
c
c  routines:
c          routine    source    brief description
c          -------------------------------------------------------------
c          i4unit     ADAS      fetch unit number for output of messages
c          i4fctn     ADAS      convert string to integer form
c          xfelem     ADAS      return element name given nuclear charge
c          xxword     ADAS      extract position of number in buffer
c          xxslen     ADAS      find string less front and tail blanks
c          xxcase     ADAS      convert a string to upper or lower case
c          xxrptn     ADAS      analyse an adf11 file partition block
c
c author:  H. P. Summers, University of Strathclyde
c          ja7.08
c          tel. +44-(0)141-548-4196
c
c date:    04/10/06
c
c version: 1.1                          date: 04/10/2006
c modified: Hugh Summers
c               - first edition.
c
c version: 1.2                          date: 21/01/2007
c modified: Allan Whiteford
c               - Commented out warning about lack of iclass,
c                 all of the present ADAS files do not contain
c                 this information
c                 (first commit to CVS)
c
c version: 1.3                          date: 08/03/2007
c modified: Hugh Summers
c               - adjustments for revised ecd formats.
c                 charge exchange donor/donor mass checks and
c                 dnr_ele, dnr_ams added to parameter return.
c
c-----------------------------------------------------------------------
      integer   nddash    , idword      , ndstack  , ndonors
c-----------------------------------------------------------------------
      character csrch1*4   , csrch2*4     , csrch3*4
      character csrch4*4   , csrch5*4
      character cdash*8    , cpart*3
c-----------------------------------------------------------------------
      parameter ( nddash = 6  , idword = 256 , ndstack = 40 )
      parameter ( csrch1 = 'iprt' , csrch2 = 'igrd'  ,csrch3 = '---/')
      parameter ( csrch4 = 'ispp' , csrch5 = 'ispb')
      parameter ( cdash = '--------' , cpart = '//#' )
      parameter ( ndonors = 4 )
c-----------------------------------------------------------------------
      integer   iunit     , i4unit   , i4fctn
      integer   iz0       , is1min   , is1max
      integer   iddimd    , itdimd   , isdimd
      integer   ndptnl    , ndptn    , ndptnc   , ndcnct
      integer   iblmx     , ismax    , itmax    , idmax
      integer   ischk
      integer   iclass    , iclass_file
      integer   i         , ic       , it         , id
      integer   icptn     , iabt     , iabt1      , iabt2    ,
     &          j         , ndash_line
      integer   nptnl     , ncnct    , ncptn_stack
      integer   nfirst    , iwords   , nwords   , ifirst     , ilast
c-----------------------------------------------------------------------
      real*8    dnr_ams   , dmass
c-----------------------------------------------------------------------
      character cstrg*80  , cstrgl*80, cterm*80   ,  chindi*4
      character xfelem*12 , str12*12 , dnr_ele*12
c-----------------------------------------------------------------------
      logical   lres      , lstan     , lptn      , lresol
      logical   lptn_old  , lwarn     , ldonor    , ldmass
c-----------------------------------------------------------------------
      integer   idash_linea(nddash)
      integer   nptn(ndptnl)          , nptnc(ndptnl,ndptn)
      integer   iptnla(ndptnl)        , iptna(ndptnl,ndptn)
      integer   iptnca(ndptnl,ndptn,ndptnc)
      integer   icnctv(ndcnct)
      integer   isstgr(isdimd)
      integer   ifirsta(idword)       , ilasta(idword)
      integer   ispbr(isdimd)         , isppr(isdimd)
c-----------------------------------------------------------------------
      real (kind=R8) ::    ddens(iddimd)         , dtev(itdimd)
      real (kind=R8) ::    drcof(isdimd,itdimd,iddimd)
c-----------------------------------------------------------------------
      character cclass(12)*4          , cpatrn(12)*4
      character cptrn1(12)*4          , cptrn2(12)*4
      character cptn_stack(ndstack)*80
      character cdonors(ndonors)*12
c-----------------------------------------------------------------------
      data      cterm /'------------------------------------------------
     &--------------------------------'/
      data      cclass/'/ACD','/SCD','/CCD','/PRB','/PRC',
     &                 '/QCD','/XCD','/PLT','/PLS','/ZCD',
     &                 '/YCD','/ECD'/
      data      cptrn1/'iprt','iprt','iprt','iprt','iprt',
     &                 'igrd','iprt','igrd','igrd','igrd',
     &                 'igrd','igrd'/
      data      cptrn2/'ispp','ispp','ispp','ispp','ispp',
     &                 'ispb','ispp','ispb','ispb','ispb',
     &                 'ispb','ispb'/
      data      cdonors/'HYDROGEN    ','DEUTERIUM   ','TRITIUM     ',
     &                  'HELIUM      '/
c-----------------------------------------------------------------------
      external  subini, subend, xerrab
c
c-----------------------------------------------------------------------
c  search for 'prnt' and count number of dash delimiters
c-----------------------------------------------------------------------

*   ..subprogram start-up calls
      call subini ('xxdata_11')

       ic = 0
       icptn = 0
       ndash_line = 0

       lstan  = .false.
       lres = .false.
       lptn   = .false.
       nptnl  = 0

   10  read(iunit,'(1a80)') cstrg
       ic = ic + 1


       call xxcase(cstrg,cstrgl,'lc')

       if((index(cstrgl,csrch1).le.0) .and.
     &    (index(cstrgl,csrch2).le.0) .and.
     &    (index(cstrgl,csrch4).le.0) .and.
     &    (index(cstrgl,csrch5).le.0) .and.
     &    (index(cstrgl,csrch3).le.0)) then

            if(index(cstrgl,cdash).gt.0) then
                ndash_line = ndash_line+1
                idash_linea(ndash_line)=ic
            endif
            if((cstrgl(1:3).eq.cpart).and.(.not.lptn)) then
                lptn = .true.
                icptn = ic
            endif
            go to 10
       else
            if((iclass.eq.1).or.(iclass.eq.3).or.
     &         (iclass.eq.4).or.(iclass.eq.5).or.
     &         (iclass.eq.6).or.(iclass.eq.7)) then
                if(ndash_line.eq.1) then
                    lstan = .true.
                elseif ((ndash_line.eq.2).and.lptn)then
                    lstan = .true.
                elseif ((ndash_line.eq.2).and.(.not.lptn))then
                    lres = .true.
                elseif ((ndash_line.eq.3).and.lptn)then
                    lres = .true.
                elseif ((ndash_line.eq.3).and.(.not.lptn))then
                    write(i4unit(-1),2002) 'header lines faulty'
                    write(i4unit(-1),2003)
                    call xerrab ('Faulty header lines')
                elseif (ndash_line.gt.3) then
                    write(i4unit(-1),2002) 'header lines faulty'
                    write(i4unit(-1),2003)
                    call xerrab ('Faulty header lines')
                endif
            elseif((iclass.eq.2).or.(iclass.eq.8).or.
     &             (iclass.eq.9).or.(iclass.eq.10).or.
     &             (iclass.eq.11).or.(iclass.eq.12)) then
                if(ndash_line.eq.1) then
                    lstan = .true.
                elseif ((ndash_line.eq.2).and.lptn)then
                    lstan = .true.
                elseif ((ndash_line.eq.2).and.(.not.lptn))then
                    lres = .true.
                elseif ((ndash_line.eq.3).and.lptn)then
                    lres = .true.
                elseif ((ndash_line.eq.3).and.(.not.lptn))then
                    write(i4unit(-1),2002) 'header lines faulty'
                    write(i4unit(-1),2003)
                    call xerrab ('Faulty header lines')
                elseif (ndash_line.gt.3) then
                    write(i4unit(-1),2002) 'header lines faulty'
                    write(i4unit(-1),2003)
                    call xerrab ('Faulty header lines')
                endif
            endif
       endif
c-----------------------------------------------------------------------
c  assign look-up names according to partitioning & set warning logical
c-----------------------------------------------------------------------

       do i=1,12
         if(lptn) then
             cpatrn(i)=cptrn2(i)
         else
             cpatrn(i)=cptrn1(i)
         endif
       enddo

       lwarn = .false.

c-----------------------------------------------------------------------
c  rewind, separate and analyse header sections according to type
c-----------------------------------------------------------------------

       rewind(iunit)

       ic =0

c-----------------------------------------------------------------------
c  first line
c-----------------------------------------------------------------------
       read(iunit,'(a80)') cstrg
       call xxcase(cstrg,cstrgl,'uc')
       ic=ic+1

       iabt = 0
       iz0 = i4fctn(cstrg(1:5),iabt)
       if(iabt.gt.0) then
            write(i4unit(-1),2002) 'faulty nuclear charge'
            write(i4unit(-1),2003)
            call xerrab ('Faulty nuclear charge')
       endif
       str12=xfelem(iz0)
       call xxslen(str12,ifirst,ilast)
       if(index(cstrgl,str12(ifirst:ilast)).le.0) then
              write(i4unit(-1),2000) 'inconsistent element name '//
     &                               'for iz0= ',iz0
              write(i4unit(-1),2001)
       endif

       iabt = 0
       idmax = i4fctn(cstrg(6:10),iabt)
       if(iabt.gt.0.or.idmax.le.0.or.idmax.gt.iddimd) then
              write(i4unit(-1),2002) 'invalid number of ','densities',
     &                                idmax,' = 0 or > ',iddimd
              write(i4unit(-1),2003)
            call xerrab ('Invalid number of densities')
       endif
       iabt = 0
       itmax = i4fctn(cstrg(11:15),iabt)
       if(iabt.gt.0.or.itmax.le.0.or.itmax.gt.itdimd) then
            write(i4unit(-1),2002) 'invalid number of ','temperatures',
     &                              itmax,' = 0 or > ',itdimd
            write(i4unit(-1),2003)
            call xerrab ('Invalid number of temperatures')
       endif
       iabt1 = 0
       is1min = i4fctn(cstrg(16:20),iabt1)
       iabt2 = 0
       is1max = i4fctn(cstrg(21:25),iabt2)
       iabt = iabt1 + iabt2

       iclass_file = 0
       do i=1,12
         if(index(cstrgl,cclass(i)).gt.0) then
            iclass_file=i
         endif
       enddo
       if(iclass_file.le.0) then
            write(i4unit(-1),2000) 'no class name given in file'
            write(i4unit(-1),2001)
       elseif(iclass_file.ne.iclass) then
            write(i4unit(-1),2004) 'file class name ',
     &                              cclass(iclass_file),
     &                             ' does not match iclass = ',iclass
            write(i4unit(-1),2003)
            call xerrab ('File class name does not match')
       endif

c-----------------------------------------------------------------------
c  check if CX donor is explicitly named in dataset top line
c-----------------------------------------------------------------------
       ldonor  = .false.
       ldmass  = .false.
       dnr_ele = ' '
       dnr_ams = 0.0_R8
       if((iclass_file.eq.3).or.(iclass_file.eq.5)) then
           do i=1,ndonors
             call xxslen(cdonors(i),ifirst,ilast)
             if(index(cstrgl,':'//cdonors(i)(ifirst:ilast)).gt.0) then
                 ldonor=.true.
                 dnr_ele=cdonors(i)
             endif
           enddo
       endif

       write(i4unit(-1),*) 'xxdata_11.for:iz0,is1min,is1max=',
     &                     iz0,is1min,is1max
c-----------------------------------------------------------------------
c  read connection vector if partial file
c-----------------------------------------------------------------------

       if(lres) then
           do i=2,idash_linea(1)
             read(iunit,'(1a80)') cstrg
           enddo

           ncnct = 0

           do i=idash_linea(1)+1,idash_linea(2)-1
             cstrg=' '
             read(iunit,'(1a80)') cstrg
             nfirst = 1
             iwords = idword
             call xxword( cstrg , ' '    , nfirst ,
     &                    iwords  ,
     &                    ifirsta , ilasta , nwords
     &                  )
             do j=1,nwords
               ncnct = ncnct + 1
               read(cstrg(ifirsta(j):ilasta(j)),*) icnctv(ncnct)
             enddo
           enddo

       else

           ncnct = 0

       endif

c
       rewind(iunit)

c-----------------------------------------------------------------------
c  read partition data if present
c-----------------------------------------------------------------------
c
       lptn_old = lptn
       if(lptn_old) then
          do i=1,icptn-1
            read(iunit,'(1a80)') cstrg
          enddo

           lresol = .false.
           call  xxrptn( iunit  , ndstack,
     &                   ndptnl , ndptn  , ndptnc ,
     &                   nptnl  , nptn   , nptnc  ,
     &                   iptnla , iptna  , iptnca ,
     &                   lresol , lptn   ,
     &                   cstrg  ,
     &                   ncptn_stack     , cptn_stack
     &                 )

       endif

c------- now check is1min and is1max

       if(.not.lptn) then

           if((iabt.gt.0).or.(is1min.gt.is1max).or.(is1min.lt.1)
     &         .or.(is1max.gt.iz0)) then
             write(i4unit(-1),2002) 'incorrect ion or partition limits'
             write(i4unit(-1),2003)
             call xerrab ('Incorrect ion or partition limits')
           endif

       elseif(lstan.and.lptn) then

           if((iabt.gt.0).or.(is1min.gt.is1max).or.(is1min.ne.1)
     &         .or.(is1max.ne.nptn(1)-1)) then
             write(i4unit(-1),2002) 'incorrect ion or partition limits'
             write(i4unit(-1),2003)
             call xerrab ('Incorrect ion or partition limits')
           endif

       elseif(lres.and.lptn) then

           if((iabt.gt.0).or.(is1min.gt.is1max).or.(is1min.ne.1)
     &         .or.(is1max.ne.ncnct-1)) then
             write(i4unit(-1),2002) 'incorrect ion or partition limits'
             write(i4unit(-1),2003)
             call xerrab ('Incorrect ion or partition limits')
           endif

       endif

c-----------------------------------------------------------------------
c  read temperatures and densities
c-----------------------------------------------------------------------

       rewind(iunit)

       do i=1,idash_linea(ndash_line)
          read(iunit,'(1a80)') cstrg
       enddo

       read(iunit,1000) ( ddens(i) , i = 1 , idmax )
       read(iunit,1000) ( dtev(i)  , i = 1 , itmax )

c------------------------------------------
c  read parent and base metastable indices
c------------------------------------------

       chindi = cpatrn(iclass)

       iblmx  = 0
   20  read(iunit,'(a80)',end=30) cstrg
       call xxcase(cstrg,cstrgl,'lc')
       if(cstrgl(2:80).ne.cterm(2:80).and.cstrgl(2:2).ne.' ') then
c           if(lptn) then
           if(lres) then
               if( cstrgl(24:27) .eq. chindi) then
                  iblmx = iblmx  + 1
                  if((iclass.eq.4).or.(iclass.eq.5).or.
     &               (iclass.eq.8).or.(iclass.eq.9).or.
     &               (iclass.eq.10).or.(iclass.eq.11))then
                      read(cstrgl,1003) isppr(iblmx),
     &                                  isstgr(iblmx)
                      ispbr(iblmx)=0
                  else
                      read(cstrgl,1002) isppr(iblmx), ispbr(iblmx),
     &                                 isstgr(iblmx)
                  endif
               elseif(( cstrgl(24:27). eq. cptrn1(iclass)).or.
     &                ( cstrgl(24:27). eq. cptrn2(iclass))) then
                  iblmx = iblmx  + 1
                  if((iclass.eq.4).or.(iclass.eq.5).or.
     &               (iclass.eq.8).or.(iclass.eq.9).or.
     &               (iclass.eq.10).or.(iclass.eq.11))then
                      read(cstrgl,1003) isppr(iblmx),
     &                                  isstgr(iblmx)
                      ispbr(iblmx)=0
                  else
                     read(cstrgl,1002) isppr(iblmx), ispbr(iblmx),
     &                                 isstgr(iblmx)
                  endif
                  if(.not.lwarn) then
                      write(i4unit(-1),2006) 'incorrect pointer code',
     &                                       ' for class: actual = ',
     &                                         cstrgl(24:27),
     &                                       ', expected = ',chindi
                      write(i4unit(-1),2001)
                      lwarn = .true.
                  endif
               else
                  write(i4unit(-1),2005) 'incorrect pointer code for',
     &                                   ' class: actual = ',
     &                                     cstrgl(24:27),
     &                                   ', expected = ',chindi
                  write(i4unit(-1),2003)
                  call xerrab ('Incorrect pointer code')
              endif
           else
              iblmx = iblmx + 1
              read(cstrgl,1001) isstgr(iblmx)
              isppr(iblmx) = 1
              if((iclass.eq.4).or.(iclass.eq.5).or.
     &           (iclass.eq.8).or.(iclass.eq.9).or.
     &           (iclass.eq.10).or.(iclass.eq.11))then
                  ispbr(iblmx) = 0
              else
                  ispbr(iblmx) = 1
              endif
           endif
c----------------------------------------------
c  get the donor mass for classes ccd and prc
c----------------------------------------------
           dmass   = 0.0_R8
           if((iclass.eq.3).or.(iclass.eq.5))then
              if((cstrgl(46:47). eq. 'mh').or.
     &           (cstrgl(46:47). eq. 'md')) then
                      read(cstrgl(49:52),'(f4.2)') dmass
              endif
              if((dmass.gt.0.0_R8).and.(dnr_ams.eq.0.0_R8))then
                  dnr_ams=dmass
                  ldmass = .true.
              elseif((dmass.gt.0.0_R8).and.(dmass.eq.dnr_ams))then
                  continue
              else
                  ldmass  = .false.
                  dnr_ams = 0.0_R8
              endif
           endif

c--------------------------
c  read final gcr values
c--------------------------

           do 25 it = 1 , itmax
             if(iclass.le.9) then
                 read(iunit,1000) ( drcof(iblmx,it,id) , id = 1 , idmax)
             else
                 read(iunit,1004) ( drcof(iblmx,it,id) , id = 1 , idmax)
             endif
   25      continue
           go to 20
       endif
c
   30  close(iunit)

c-----------------------------------------------------------------------
c  issue warnings if donor element and/or donor mass is ambiguous
c-----------------------------------------------------------------------
       if(((iclass.eq.3).or.(iclass.eq.5)).and.(.not.ldonor)) then
           write(i4unit(-1),2000) 'unspecified CX donor element'
           write(i4unit(-1),2001)
       endif
       if(((iclass.eq.3).or.(iclass.eq.5)).and.(.not.ldmass)) then
           write(i4unit(-1),2000) 'unspecified CX donor mass'
           write(i4unit(-1),2001)
       endif
c
c-----------------------------------------------------------------------
c      verify s1 set in master file consistent with is1min and is1max
c      except for xcd and qcd cases
c-----------------------------------------------------------------------
c
       if(iclass.eq.6.or.iclass.eq.7) then
           ismax = is1max-is1min+1
       else
           if(iclass.eq.12) then
               ischk=is1min-1
           else
               ischk = is1min
           endif

           do i=1,iblmx
             if(isstgr(i).ne.ischk) then
                 ischk = ischk+1
             endif
          enddo

          if(((iclass.ne.10).and.(iclass.ne.11).and.(iclass.ne.12))
     &        .and.(ischk.eq.is1max)) then
               ismax = is1max-is1min+1
          elseif(((iclass.eq.10).or.(iclass.eq.11)).and.
     &        (ischk.eq.is1max+1)) then
               ismax = is1max-is1min+2
          elseif((iclass.eq.12).and.
     &        (ischk.eq.is1max)) then
               ismax = is1max-is1min+1
               ismax = is1max-is1min+2
          else
                write(i4unit(-1),2002) 'inconsistent s1 set in file'
                write(i4unit(-1),2003)
                call xerrab ('Inconsistent s1 set')
           endif
       endif
c
c-----------------------------------------------------------------------
c     make up connection vectors and root partitions if appropriate
c     in the resolved case with 0 root partition level, the sum over
c     the connection vector gives the number of partitions
c-----------------------------------------------------------------------
c
      if ((.not.lres).and.(.not.lptn)) then
          ncnct = iz0+1
          do i=1,ncnct
             icnctv(i)=1
          enddo
      elseif ((.not.lres).and.lptn) then
          ncnct = is1max+1
          do i=1,ncnct
             icnctv(i)=1
          enddo
      endif

      if ((.not.lres).and.(.not.lptn)) then
          nptnl   = 1
          iptnla(nptnl) = 1
          nptn(1) =  iz0+1
          do i=1,nptn(1)
            iptna(nptnl,i) = i-1
            nptnc(nptnl,i) = 1
            iptnca(1,iptna(nptnl,i)+1,nptnc(nptnl,i))=i-1
          enddo
      elseif (lres.and.(.not.lptn)) then

          nptnl   = 1
          iptnla(nptnl) = 0
          nptn(1)=1
          do i=1,ncnct
            nptn(1)=nptn(1)+icnctv(i)
          enddo
          do i=1,nptn(1)
            iptna(nptnl,i) = i-1
            nptnc(nptnl,i) = 1
            iptnca(1,iptna(nptnl,i)+1,nptnc(nptnl,i))=i-1
          enddo
      endif
c
c-----------------------------------------------------------------------
c
 1000 format(8f10.5)
 1001 format(57x,i2)
 1002 format(28x,i2,9x,i2,16x,i2)
 1003 format(28x,i2,9x,2x,16x,i2)
 1004 format(1p8d10.3)
c
 2000 format(1x,30('*'),' xxdata_11 warning ',30('*')//
     &       2x,a,a,i3,a,i3 )
 2001 format(/1x,30('*'),' program continues ',30('*'))
 2002 format(1x,30('*'),'  xxdata_11 error  ',30('*')//
     &       2x,a,a,i3,a,i3 )
 2003 format(/1x,30('*'),' program terminated ',29('*'))
 2004 format(1x,30('*'),'  xxdata_11 error  ',30('*')//
     &       2x,a,a,a,i3)
 2005 format(1x,30('*'),'  xxdata_11 error  ',30('*')//
     &       2x,a,a,a4,a,a4 )
 2006 format(1x,30('*'),' xxdata_11 warning ',30('*')//
     &       2x,a,a,a4,a,a4 )
c
c-----------------------------------------------------------------------
c
      call subend ()
      return
      end subroutine xxdata_11
CX ULTRIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adaslib/utility/i4fctn.for,v 1.3 2004/07/06 14:06:47 whitefor Exp $ Date $Date: 2004/07/06 14:06:47 $
CX
      FUNCTION I4FCTN( STR , IABT )
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  *************** FORTRAN77 INTEGER*4 FUNCTION: I4FCTN ****************
C
C  FUNCTION: TO CONVERT AN INTEGER NUMBER STORED IN THE STRING 'STR'
C            INTO A INTEGER*4 VARIABLE, USING INTERNAL READ.
C            INTIALLY THE PROGRAM CHECKS TO SEE IF THE NUMBER IS OF A
C            VALID FORM.
C
C  CALLING PROGRAM: GENERAL USE
C
C  FUNCTION:
C
C          (I*4)   I4FCTN  = FUNCTION NAME
C          (C*(*)) STR     = STRING CONTAINING SINGLE FLOATING POINT NO.
C          (I*4)   IABT    = RETURN CODE:
C                               0 => NO ERROR
C                               1 => ERROR (A VALUE 'I4FCTN=0' WILL BE
C                                           RETURNED).
C
C          (C*1)   CH0      = PARAMETER = '0'
C          (C*1)   CH9      = PARAMETER = '9'
C          (C*1)   BLANK    = PARAMETER = ' '
C          (C*1)   CPLUS    = PARAMETER = '+'
C          (C*1)   CMINUS   = PARAMETER = '-'
C
C          (I*4)   ILEN     = LENGTH OF 'STR' STRING IN BYTES
C          (I*4)   ILAST    = POSITION OF LAST BYTE OF IDENTIFIED NUMBER
C          (I*4)   I1       = STARTING BYTE IN 'STR' OF NUMBER
C                             INCLUDING SIGN IF PRESENT
C          (I*4)   IS       = 0 => NUMBER HAS NO SIGN
C                             1 => NUMBER HAS A SIGN
C          (I*4)   ICH0     = ICHAR('0')
C          (I*4)   ICH9     = ICHAR('9')
C          (I*4)   ISTR     = ICHAR(CURRENT BYTE POSITION IN 'STR')
C          (I*4)   I        = GENERAL USE
C
C          (L*4)   LFOUND   = .TRUE.  => ALL OF THE INPUT NUMBER BYTES
C                                        HAVE BEEN ASSESSED.
C                             .FALSE. => INPUT NUMBER BYTES STILL BEING
C                                        ASSESSED.
C          (L*4)   LSTART   = .TRUE.  => THE FIRST DIGIT HAS BEEN FOUND
C                             .FALSE. => THE FIRST DIGIT HAS NOT YET
C                                        BEEN REACHED.
C
C          (C*5)   CFORM5   = FORMAT FOR INTERNAL READING OF INTEGER
C
C
C NOTE:     AN ERROR WILL OCCUR (IABT=1) IF THERE IS MORE THAN ONE
C           NUMBER OCCURING IN THE STRING 'STR()'
C
C
C AUTHOR:   PAUL E. BRIDEN (TESSELLA SUPPORT SERVICES PLC)
C           K1/0/37
C           JET EXT. 2520
C
C DATE:     11/07/90
C
C UPDATE:   11/02/92 - PE BRIDEN: BLANKS NOW ALLOWED BETWEEN SIGN AND
C                                 FIRST DIGIT. LSTART VARIABLE ADDED.
C                                 VARIABLE I2 REMOVED.
C                                 + SOME MINOR RECODING - (IF STRING
C                                 ENTERED IS BLANK IABT IS NOW SET TO 1)
C
C UPDATE:   16/08/93 - PE BRIDEN: CORRECTED BUG TO ALLOW BLANKS BETWEEN
C                                 SIGN AND FIRST DIGIT (SEE ABOVE).
C                                 1) ILAST VARIABLE ADDED.
C                                 2) FORMATTED READ USED INSTEAD OF *
C                                    WHEN CONVERTING IDENTIFIED INTEGER
C                                    USING THE INTERNAL READ. (THIS
C                                    RESTRICTS IDENTIFIED NUMBER TO BE
C                                    < 100 BYTES IN LENGTH!)
C                                 3) EXCLUDE TRAILING BLANKS IN THE
C                                    INTERNAL READING OF THE INTEGER
C                                    I.E. STR(I1:ILAST) INSTEAD OF
C                                         STR(I1:ILEN)
C
C UPDATE:   07/03/95 - PE BRIDEN: INSTEAD OF USING FORMAT SPECIFIER I99
C                                 WHEN INTERNALLY READING THE INTEGER
C                                 CREATE THE APPROPRIATE SPECIFIER
C                                 WITHIN CFORM5 AND USE THIS.
C
C VERSION  : 1.3
C DATE     : 20-12-2001
C MODIFIED : Martin O'Mullane
C               - Removed mainframe listing information beyond column 72.
C
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      CHARACTER  CH0*1   , CH9*1   , BLANK*1   , CPLUS*1   , CMINUS*1
C-----------------------------------------------------------------------
      PARAMETER( CH0='0', CH9='9', BLANK=' ', CPLUS='+', CMINUS='-' )
C-----------------------------------------------------------------------
      CHARACTER  STR*(*)
C-----------------------------------------------------------------------
      INTEGER    I4FCTN  , IABT
      INTEGER    I1      , IS      , ILEN      , ILAST     ,
     &           ICH0    , ICH9    , ISTR      , I
C-----------------------------------------------------------------------
      LOGICAL    LSTART  , LFOUND
C-----------------------------------------------------------------------
      CHARACTER  CFORM5*5
C-----------------------------------------------------------------------
      DATA       CFORM5 / '(I??)' /
C-----------------------------------------------------------------------
      SAVE       CFORM5
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C INITIALISE VALUES
C-----------------------------------------------------------------------
C
      I4FCTN = 0
      IABT   = 0
      I1     = 0
      IS     = 0
      LSTART = .FALSE.
      LFOUND = .FALSE.
      ICH0   = ICHAR(CH0)
      ICH9   = ICHAR(CH9)
      ILEN   = LEN(STR)
      ILAST  = ILEN
C
C-----------------------------------------------------------------------
C FIND STARTING BYTE OF NUMBER
C-----------------------------------------------------------------------
C
         DO 1 I=1,ILEN
               IF ( STR(I:I).NE.BLANK ) THEN
                  I1 = I
                  GOTO 2
               ENDIF
    1    CONTINUE
C
C-----------------------------------------------------------------------
C IDENTIFY IF NUMBER HAS A SIGN
C-----------------------------------------------------------------------
C
    2    IF (I1.EQ.0) THEN
            IABT = 1
            RETURN
         ENDIF
C
      IF   ( ( STR(I1:I1).EQ.CPLUS  )
     &                  .OR.
     &       ( STR(I1:I1).EQ.CMINUS ) ) IS=1
C
C-----------------------------------------------------------------------
C IDENTIFY IF NUMBER IS OF A VALID FORM
C-----------------------------------------------------------------------
C
         DO 3 I=I1+IS,ILEN

               IF (LFOUND) THEN
C
C-----------------------------------------------------------------------
C INPUT NO. COMPLETELY DEFINED: IDENTIFY IF EXTRA NON-BLANK BYTES EXIST
C-----------------------------------------------------------------------
C
                  IF (STR(I:I).NE.BLANK) IABT=1
C-----------------------------------------------------------------------
               ELSEIF (STR(I:I).EQ.BLANK) THEN
                  LFOUND = LSTART
C-----------------------------------------------------------------------
               ELSE
                  LSTART = .TRUE.
                  ILAST  = I
                  ISTR   = ICHAR(STR(I:I))
                  IF ( (ISTR.LT.ICH0) .OR. (ISTR.GT.ICH9) ) IABT=1
               ENDIF
C
C-----------------------------------------------------------------------
C RETURN ERROR CODE IF ERROR FOUND
C-----------------------------------------------------------------------
C
            IF (IABT.NE.0) RETURN
    3    CONTINUE
C
C-----------------------------------------------------------------------
C IDENTIFY IF VALID NUMBER FOUND (RECODED: PEB 11/02/92)
C                                (RECODED: PEB 07/03/95 - ADDED CFORM5)
C YES => USE INTERNAL READ TO OBTAIN THE INTEGER NUMBER
C NO  => RETURN ERROR CODE IF ERROR FOUND
C-----------------------------------------------------------------------
C
         IF (LSTART) THEN
            I      = 1 + ILAST - I1
            I      = MIN0(I,99)
            WRITE(CFORM5(3:4),'(I2.2)') I
            READ(STR(I1:ILAST),CFORM5) I4FCTN
         ELSE
            IABT=1
         ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END FUNCTION I4FCTN
CX UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adaslib/system/i4unit.for,v 1.1 2004/07/06 14:08:05 whitefor Exp $ Date $Date: 2004/07/06 14:08:05 $
CX
      FUNCTION I4UNIT( IUNIT )
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  ************** FORTRAN77 INTEGER*4 FUNCTION: I4UNIT *****************
C
C  PURPOSE: TO RESET OR RETURN A STORED INTEGER*4 VALUE GREATER THAN OR
C           EQUAL TO ZERO.
C           THIS IS USED WITHIN ADAS TO STORE THE STREAM/UNIT NUMBER
C           FOR THE OUTPUT OF ERROR MESSAGES (TO THE SCREEN).
C
C           BY DEFAULT THE STORED VALUE WILL BE 6, AND WILL BE RETURNED
C           BY THE FUNCTION IF IUNIT ON INPUT < 0.
C
C           TO RESET THE STORED VALUE THEN SET IUNIT TO THE REQUIRED
C           POSITIVE INTEGER (INC. ZERO). THIS VALUE WILL ALSO BE
C           RETURNED BY THE FUNCTION.
C
C                 IUNIT VALUE               RETURNED FUNCTION VALUE
C                 -----------               -----------------------
C                 IUNIT <  0            = CURRENT STORED INTEGER VALUE
C                                         (6 BY DEFAULT).
C                 IUNIT >= 0            = IUNIT , AND RESETS THE STORED
C                                                 VALUE TO IUNIT.
C
C
C  CALLING PROGRAM: GENERAL USE
C
C  SUBROUTINE:
C
C  O     : (I*4)  I4UNIT   = FUNCTION NAME - (SEE ABOVE)
C
C  I     : (I*4)  IUNIT    = FUNCTION ARGUMENT - (SEE ABOVE)
C
C          (I*4)  IDEFLT   = PARAMETER = DEFAULT STORED INTEGER VALUE
C
C          (I*4)  ICURNT   = CURRENT STORED INTEGER VALUE
C
C
C ROUTINES:
C          ROUTINE    SOURCE    BRIEF DESCRIPTION
C          ------------------------------------------------------------
C
C
C AUTHOR:  PAUL E. BRIDEN (TESSELLA SUPPORT SERVICES PLC)
C          K1/0/37
C          JET EXT. 5023
C
C DATE:    23/04/93
C
C UPDATE:  24/05/93 - PE BRIDEN - ALLOWED 0 TO BE A VALID STORED NUMBER
C
C-----------------------------------------------------------------------
      INTEGER    I4UNIT     , IDEFLT
C-----------------------------------------------------------------------
      PARAMETER( IDEFLT = 6 )
C-----------------------------------------------------------------------
      INTEGER    IUNIT      , ICURNT
C-----------------------------------------------------------------------
      SAVE       ICURNT
C-----------------------------------------------------------------------
      DATA       ICURNT / IDEFLT /
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C RETRIEVE OR RESET STORED INTEGER VALUE ACCORDINGLY
C-----------------------------------------------------------------------
C
      IF (IUNIT.LT.0) THEN
        I4UNIT = ICURNT
      ELSE
        ICURNT = IUNIT
        I4UNIT = IUNIT
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END FUNCTION I4UNIT
CX UNIX PORT - SCCS Info : Module @(#)$Header: /home/adascvs/fortran/adaslib/atomic/xfelem.for,v 1.2 2004/07/06 15:30:02 whitefor Exp $ Date $Date: 2004/07/06 15:30:02 $
CX
      FUNCTION XFELEM ( IZ0 )
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C ************* FORTRAN77 CHARACTER*12 FUNCTION: XFELEM ****************
C
C PURPOSE: TO RETURN THE NAME OF THE ELEMENT WITH NUCLEAR CHARGE IZ0
C          (CHARACTER*12 FUNCTION VERSION OF 'XXELEM')
C
C CALLING PROGRAM: GENERAL USE
C
C FUNCTION:
C
C          (C*12) XFELEM  = FUNCTION NAME -
C                           NAME OF ELEMENT WITH NUCLEAR CHARGE 'IZ0'
C          (I*4)  IZ0     = ELEMENT NUCLEAR CHARGE
C
C          (C*12) NAMES() = NAMES OF FIRST 50 ELEMENTS.
C                           ARRAY DIMENSION => NUCLEAR CHARGE
C
C NOTES:    IF NUCLEAR CHARGE IS OUT OF RANGE, I.E.NOT BETWEEN 1 & 50,
C           THEN THE CHARACTER STRING 'XFELEM' IS RETURNED BLANK.
C
C ROUTINES: NONE
C
C
C AUTHOR:   PAUL E. BRIDEN (TESSELLA SUPPORT SERVICES PLC)
C           K1/0/81
C           JET EXT. 4569
C
C DATE:     13/02/91
C
C VERSION: 1.2
C UPDATE:   17/09/99  HUGH SUMMERS - INCREASED ELEMENT NUMBER TO 92
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      INTEGER      IZ0
C-----------------------------------------------------------------------
      CHARACTER*12 XFELEM , NAMES(92)
C-----------------------------------------------------------------------
      DATA NAMES/'HYDROGEN    ','HELIUM      ','LITHIUM     ',
     &           'BERYLLIUM   ','BORON       ','CARBON      ',
     &           'NITROGEN    ','OXYGEN      ','FLUORINE    ',
     &           'NEON        ','SODIUM      ','MAGNESIUM   ',
     &           'ALUMINIUM   ','SILICON     ','PHOSPHORUS  ',
     &           'SULPHUR     ','CHLORINE    ','ARGON       ',
     &           'POTASSIUM   ','CALCIUM     ','SCANDIUM    ',
     &           'TITANIUM    ','VANADIUM    ','CHROMIUM    ',
     &           'MANGANESE   ','IRON        ','COBALT      ',
     &           'NICKEL      ','COPPER      ','ZINC        ',
     &           'GALLIUM     ','GERMANIUM   ','ARSENIC     ',
     &           'SELENIUM    ','BROMINE     ','KRYPTON     ',
     &           'RUBIDIUM    ','STRONTIUM   ','YTTRIUM     ',
     &           'ZIRCONIUM   ','NIOBIUM     ','MOLYBDENUM  ',
     &           'TECHNETIUM  ','RUTHENIUM   ','RHODIUM     ',
     &           'PALLADIUM   ','SILVER      ','CADMIUM     ',
     &           'INDIUM      ','TIN         ','ANTIMONY    ',
     &           'TELLURIUM   ','IODINE      ','XENON       ',
     &           'CESIUM      ','BARIUM      ','LANTHANUM   ',
     &           'CERIUM      ','PRAESODYMIUM','NEODYMIUM   ',
     &           'PROMETHIUM  ','SAMARIUM    ','EUROPIUM    ',
     &           'GADOLINIUM  ','TERBIUM     ','DYSPROSIUM  ',
     &           'HOLMIUM     ','ERBIUM      ','THULIUM     ',
     &           'YTTERBIUM   ','LUTETIUM    ','HAFNIUM     ',
     &           'TANTALUM    ','TUNGSTEN    ','RHENIUM     ',
     &           'OSMIUM      ','IRIDIUM     ','PLATINUM    ',
     &           'GOLD        ','MERCURY     ','THALLIUM    ',
     &           'LEAD        ','BISMUTH     ','POLONIUM    ',
     &           'ASTATINE    ','RADON       ','FRANCIUM    ',
     &           'RADIUM      ','ACTINIUM    ','THORIUM     ',
     &           'PROTACTINIUM','URANIUM     '/
C-----------------------------------------------------------------------
         IF ( (IZ0.GT.92).OR.(IZ0.LT.0) ) THEN
            XFELEM = ' '
         ELSE
            XFELEM = NAMES(IZ0)
         ENDIF
C-----------------------------------------------------------------------
      RETURN
      END FUNCTION XFELEM
      subroutine xxcase(input,output,type)

      IMPLICIT NONE

C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 SUBROUTINE: XXCASE *********************
C
C  PURPOSE: Change a string of arbitrary size into all upper case
C           or all lower case
C
C  CALLING PROGRAM: GENERAL USE.
C
C  INPUT    : (C*(*)) INPUT = Input String
C  INPUT    : (C*2)   TYPE = Type of case to convert to:
C                       'UC' -> Convert to Upper Case
C                       'LC' -> Convert to Lower Case
C                       Anything else -> No conversion
C
C  OUTPUT   : (C*(*)) OUTPUT = Output string in selected case
C
C  ROUTINES : NONE
C
C  AUTHOR   : Allan Whiteford,
C             University of Strathclyde
C
C  VERSION  : 1.1
C  DATE     : 05/09/2001
C  MODIFIED : Allan Whiteford
C             First version.
C
C  VERSION  : 1.2
C  DATE     : 05/05/2005
C  MODIFIED : Martin O'Mullane
C             The routine converted length-1 rather than the whole
C             input string.
C
C-----------------------------------------------------------------------
       integer i
       integer size
C----------------------------------------------------------------------
       character*(*) input
       character*(*) output
       character*2 type
C----------------------------------------------------------------------

       size=len(input)

       write(output(1:size),'(A)') input(1:size)
       i=1

       if (type.eq.'UC' .or. type.eq.'uc') then
   10     if    ( ichar(input(i:i)) .ge. ichar('a')
     &    .and.   ichar(input(i:i)) .le. ichar('z')
     &          ) output(i:i)=char(ichar(input(i:i))+ichar('A')
     &                                              -ichar('a'))

          i=i+1
          if (i .le. size) goto 10
       endif

       if (type.eq.'LC' .or. type.eq.'lc') then
   20     if    ( ichar(input(i:i)) .ge. ichar('A')
     &    .and.   ichar(input(i:i)) .le. ichar('Z')
     &          ) output(i:i)=char(ichar(input(i:i))+ichar('a')
     &                                              -ichar('A'))

          i=i+1
          if (i .le. size) goto 20
       endif

       return
       end subroutine xxcase
       subroutine xxrptn( iunit  , ndstack,
     &                    ndptnl , ndptn  , ndptnc ,
     &                    nptnl  , nptn   , nptnc  ,
     &                    iptnla , iptna  , iptnca ,
     &                    lresol , lptn   ,
     &                    cstrg  ,
     &                    ncptn_stack     , cptn_stack
     &                  )
       implicit none
c-----------------------------------------------------------------------
c
c  ****************** fortran77 subroutine: xxprtn *********************
c
c  purpose:  to read and analyse a partition block in a datafile header
c
c  calling program: adas416
c
c  notes:  (1) partition levels, partitions and partition components are
c              labelled starting at 0 (but see (2))
c          (2) partition level 0 labels the resolved root partition level
c              partition level 1 labels the unresolved root partition
c              level.
c          (3) for an unresolved (standard) file, the partitions are each
c              ionisation stage from the neutral to the bare nucleus and
c              they are labelled by the ion charge.  Each partition has
c              just the one component
c          (4) distinguish the indexing (starting at 1) from the label
c              (starting at 0)
c
c  subroutine:
c
c  input : (i*4)  iunit     = unit to which input file is allocated
c  input : (i*4)  ndstack   = maximum no. of text lines in partition block
c
c  input : (i*4)  ndptnl    = maximum level of partitions
c  input : (i*4)  ndptn     = maximum no. of partitions in one level
c  input : (i*4)  ndptnc    = maximum no. of components in a partition
c  input : (l*4)  lresol    = .true.  => resolved root partition
c                           = .false. => standard root partition
c
c  output: (i*4)  nptnl     = number of partition levels in block
c  output: (i*4)  nptn()    = number of partitions in partition level
c                             1st dim: partition level
c  output: (i*4)  nptnc(,)  = number of components in partition
c                             1st dim: partition level
c                             2nd dim: member partition in partition level
c  output: (i*4)  iptnla()  = partition level label (0=resolved root,1=
c                                                      unresolved root)
c                             1st dim: partition level index
c  output: (i*4)  iptna(,)  = partition member label (labelling starts at 0)
c                             1st dim: partition level index
c                             2nd dim: member partition index in partition
c                             level
c  output: (i*4)  iptnca(,,)= component label (labelling starts at 0)
c                             1st dim: partition level index
c                             2nd dim: member partition index in partition
c                             level
c                             3rd dim: component index of member partition
c  output: (l*4)  lptn      = .true.  => partition block present
c                           = .false. => partition block not present
c  output: (c*80) cstrg     = string marking end of partition block
c  output: (i*4)  ncptn_stack= number of text lines in partition block
c  output: (c*80) cptn_stack()=text lines of partition block
c                              1st dim: text line pointer
c
c
c  routines:
c          routine    source    brief description
c          -------------------------------------------------------------
c          i4unit     ADAS      fetch unit number for output of messages
c          xxslen     ADAS      find non-blank characters in string
c          xxword     ADAS      extract position of number in buffer
c
c author:  H. P. Summers, University of Strathclyde
c          ja7.08
c          tel. +44-(0)141-548-4196
c
c date:    25/08/05
c
c version: 1.1                          date: 25/08/2005
c modified: Hugh Summers
c               - first edition.
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      integer   idptnl   , idword
c-----------------------------------------------------------------------
      parameter (idptnl = 4  , idword = 256 )
c-----------------------------------------------------------------------
      integer   iunit     , i4unit
      integer   ndptnl    , ndptn      , ndptnc    , ndstack
      integer   nptnl
      integer   ifirst    , ilast
      integer   nchari    , ncharf
      integer   nfirst    , iwords     , nwords
      integer   nfirstc   , iwordsc    , nwordsc
      integer   i         , j          , k         , ic
      integer   ncptn_stack
c-----------------------------------------------------------------------
      character cstrg*80
      character strg*256
c-----------------------------------------------------------------------
      logical   lresol    , lptn
      logical   lptn_in
c-----------------------------------------------------------------------
      integer   nptn(ndptnl)          , nptnc(ndptnl,ndptn)
      integer   iptnla(ndptnl)        , iptna(ndptnl,ndptn)
      integer   iptnca(ndptnl,ndptn,ndptnc)
      integer   ifirsta(idword)       , ilasta(idword)
      integer   ifirstac(idword)      , ilastac(idword)
      integer   ncpta(idptnl)
c-----------------------------------------------------------------------
      character cptna(idptnl)*1024
      character cptn_stack(ndstack)*80
c-----------------------------------------------------------------------
      external  subini, subend, xerrab
c
c-----------------------------------------------------------------------
c  read partition data if present
c-----------------------------------------------------------------------

       call subini ('xxrptn')
       if(idptnl.ne.ndptnl) then
         write(i4unit(-1),1001) 'mismatch of internal dimensions'
         write(i4unit(-1),1002)
         call xerrab ('Internal dimensions mismatch')
       endif

       do i=1,ndptnl
         cptna(i) = ' '
       enddo
c
       lptn   = .false.
       cstrg  = ' '
       nptnl  = 0
       ncptn_stack = 0

   10  read(iunit,'(1a80)') cstrg

       if(cstrg(1:3).eq.'//#') then
           nchari = 0
           nptnl = nptnl + 1
           lptn = .true.
           lptn_in  = .true.
           call xxslen(cstrg,ifirst,ilast)

           ncptn_stack=ncptn_stack+1
           cptn_stack(ncptn_stack)=' '
           cptn_stack(ncptn_stack)(1:ilast-ifirst+1)=cstrg(ifirst:ilast)

           ncharf = nchari+ilast-ifirst+1
           cptna(nptnl)(nchari+1:ncharf) = cstrg(ifirst:ilast)
           nchari=ncharf
       elseif((cstrg(1:3).eq.'   ').and.lptn_in) then
           call xxslen(cstrg,ifirst,ilast)

           ncptn_stack=ncptn_stack+1
           cptn_stack(ncptn_stack)=' '
           cptn_stack(ncptn_stack)(1:ilast-ifirst+7)=
     &                     '      '//cstrg(ifirst:ilast)
           ncharf=nchari+ilast-ifirst+2
           cptna(nptnl)(nchari+1: ncharf) = ' '//cstrg(ifirst:ilast)
           nchari = ncharf
       elseif(cstrg(1:3).eq.'---') then
           lptn_in = .false.
           go to 15
       else
           write(i4unit(-1),1001)
           write(i4unit(-1),1002)
           call xerrab ('Faulty string')
       endif
       go to 10
   15  continue

c
c-----------------------------------------------------------------------
c  preliminary checks
c-----------------------------------------------------------------------
c
       if(lptn) then
           do i = 1,nptnl
             read(cptna(i)(4:5),'(i2)') iptnla(i)
c             write(i4unit(-1),*) 'xxrptn: i,iptnla=',i,iptnla(i)
           enddo
           if((iptnla(1).gt.2).and.(nptnl.eq.1).and.(.not.lresol)) then
               write(i4unit(-1),1001) 'partition level',i,' in error'
               write(i4unit(-1),1002)
               call xerrab ('Error in partition level')
           endif
           do i = 1,nptnl
             call xxslen(cptna(i),ifirst,ilast)
             if(cptna(i)(ilast:ilast).ne.'/') then
                 write(i4unit(-1),1001) 'partition level',i,
     &                                  ' not terminated'
                 write(i4unit(-1),1002)
                 call xerrab ('Partition level not terminated')
             endif
           enddo

        endif
c
c-----------------------------------------------------------------------
c  analyse partitions
c-----------------------------------------------------------------------
c
       do i=1,nptnl
           nfirst = 1
           iwords = idword
           call xxslen(cptna(i),ifirst,ilast)
           call xxword( cptna(i)(ifirst+5:ilast) , '/'    , nfirst ,
     &                   iwords  ,
     &                   ifirsta , ilasta , nwords
     &                 )
           if((nwords-2*int(nwords/2)).ne.0) then
                 write(i4unit(-1),1001) 'partition level',i,
     &                                  ' partition count wrong'
                 write(i4unit(-1),1002)
                 call xerrab ('Wrong partition count')
           endif

           nptn(i)=nwords/2
           ic = 0
c           write(i4unit(-1),*) 'xxrptn: i,nwords,nptn(i)=',
c     &                                 i,nwords,nptn(i)
           do j=1,nptn(i)
              strg=' '
              strg=cptna(i)(ifirsta(2*j-1)+5:ilasta(2*j-1)+5)
c              write(i4unit(-1),*) 'xxrptn-xxslen:strg=',strg
              call xxslen(strg,ifirst,ilast)
              read(strg(ifirst+1:ilast),'(i2)') iptna(i,j)
              strg = ' '
              strg=cptna(i)(ifirsta(2*j)+5:ilasta(2*j)+5)

              iwordsc = idword
              nfirstc = 1
c              write(i4unit(-1),*) 'xxrptn-xxword:strg=',strg
c              write(i4unit(-1),*) 'xxrptn-xxword:nfirstc=',nfirstc
              call xxword( strg , ' '    , nfirstc ,
     &                     iwordsc  ,
     &                     ifirstac , ilastac , nwordsc
     &                   )
c              write(i4unit(-1),*) 'xxrptn:nwordsc=',nwordsc
              nptnc(i,j)=nwordsc
c              write(i4unit(-1),*) 'xxrptn: i,j,nptnc(i,j)=',
c     &                                     i,j,nptnc(i,j)
c              write(i4unit(-1),*) 'xxrptn: ifirstac(1),ilastac(1)=',
c     &                                     ifirstac(1),ilastac(1)
c            write(i4unit(-1),*) 'xxrptn:',strg(ifirstac(1):ilastac(1))

              do k=1,nptnc(i,j)
                read(strg(ifirstac(k):ilastac(k)),*) iptnca(i,j,k)
c------------------------------------------------
c  components should be sequential and contiguous.
c-------------------------------------------------
                if(iptnca(i,j,k).ne.ic) then
                     write(i4unit(-1),1001) 'partition count',ic,
     &                                      ' out of sequence'
                     write(i4unit(-1),1002)
                     call xerrab ('Partition count out of sequence')
                else
                     ic=ic+1
                endif

              enddo
           enddo

           ncpta(i)=ic
c          write(i4unit(-1),*) 'i,ncpta(i)=',i,ncpta(i)

        enddo
c
c-----------------------------------------------------------------------
c  final checks. level i+1 partition components should span
c  partition i p-values
c-----------------------------------------------------------------------
c
       if (nptnl.ge.2)then
           do i=2,nptnl
             if(ncpta(i-1).ne.nptn(i))then
                 write(i4unit(-1),1001) 'partition level',i,
     &                                  ' count incorrect'
                 write(i4unit(-1),1002)
                 call xerrab ('Incorrect partition count')
             endif
           enddo
       endif

       call subend ()
       return
c
c-----------------------------------------------------------------------
c
 1001  format(1x,32('*'),' xxrptn error ',32('*')//
     &        1x,'fault in input data file: ',a,i3,a)
 1002  format(/1x,29('*'),' program terminated ',29('*'))

       end subroutine xxrptn
CX UNIX PORT - SCCS Info : Module @(#)$Header: /home/adascvs/fortran/adaslib/utility/xxslen.for,v 1.1 2004/07/06 15:39:12 whitefor Exp $ Date $Date: 2004/07/06 15:39:12 $
CX
      SUBROUTINE XXSLEN( CSTRNG , IFIRST , ILAST )
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 SUBROUTINE: XXSLEN *********************
C
C  PURPOSE: TO IDENTIFY THE FIRST AND LAST NON-BLANK CHARACTER IN A
C           STRING. (IF INPUT STRING IS BLANK IFIRST=ILAST=0)
C
C  CALLING PROGRAM: GENERAL USE
C
C  SUBROUTINE:
C
C  INPUT : (C*(*)) CSTRNG   = INPUT STRING FOR INTERROGATION
C
C  OUTPUT: (I*4)   IFIRST   = BYTE POSITION OF FIRST NON-BLANK CHARACTER
C                             IN INPUT STRING.
C  OUTPUT: (I*4)   ILAST    = BYTE POSITION OF LAST  NON-BLANK CHARACTER
C                             IN INPUT STRING.
C
C          (I*4)   I        = GENERAL USE
C          (I*4)   ILEN     = LENGTH OF 'CSTRNG' STRING IN BYTES
C
C ROUTINES: NONE
C
C NOTE:
C
C
C AUTHOR:  PAUL E. BRIDEN (TESSELLA SUPPORT SERVICES PLC)
C          K1/0/37
C          JET EXT. 6023
C
C DATE  :  06/07/93
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      INTEGER     IFIRST     , ILAST    , ILEN    , I
C-----------------------------------------------------------------------
      CHARACTER   CSTRNG*(*)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C
      ILEN   = LEN(CSTRNG)
C-----------------------------------------------------------------------
      IFIRST = 0
      ILAST  = 0
C-----------------------------------------------------------------------
C
         DO 1 I=1,ILEN
C
            IF (CSTRNG(I:I).NE.' ') THEN
               IF (IFIRST.EQ.0) IFIRST = I
               ILAST = I
            ENDIF
C
    1    CONTINUE
C
C-----------------------------------------------------------------------
C
       RETURN
       END SUBROUTINE XXSLEN
CX UNIX PORT - SCCS Info : Module @(#)$Header: /home/adascvs/fortran/adaslib/utility/xxword.for,v 1.1 2004/07/06 15:40:43 whitefor Exp $ Date $Date: 2004/07/06 15:40:43 $
CX
      SUBROUTINE XXWORD( CTEXT  , CDELIM , NFIRST ,
     &                   IWORDS ,
     &                   IFIRST , ILAST  , NWORDS
     &                 )
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 SUBROUTINE: XXWORD *********************
C
C  PURPOSE: TO EXTRACT THE Nfirst to (Nfirst+IWORDS-1) WORDS FROM AN
C           INPUT STRING. OUTPUTS THE FIRST AND LAST BYTE INDICES OF
C           EACH WORD AS WELL AS THE TOTAL NUMBER OF WORDS FOUND.
C
C           A WORD = A STRING OF CHARACTERS SEPARATED BY ANY CHARACTER
C                    CONTAINED IN THE INPUT STRING CDELIM.
C
C  CALLING PROGRAM: GENERAL USE
C
C  SUBROUTINE:
C
C  INPUT : (C*(*)) CTEXT   = INPUT TEXT LINE CONTAINING STRING
C  INPUT : (C*(*)) CDELIM  = INPUT STRING CONTAINING DELIMITER CHARS.
C  INPUT : (I*4)   NFIRST  = THE INDEX NO. OF THE FIRST WORD TO EXTRACT.
C
C  I/O   : (I*4)   IWORDS  = INPUT : SIZE OF IFIRST, ILAST(ARRAYS)
C                                    (I.E. NUMBER OF WORDS TO EXTRACT)
C                          = OUTPUT: NUMBER OF REQUESTED WORDS FOUND
C
C  OUTPUT: (I*4)   IFIRST()= INDEX OF FIRST BYTE OF THE Nth WORD
C  OUTPUT: (I*4)   ILAST() = INDEX OF LAST  BYTE OF THE Nth WORD
C  OUTPUT: (I*4)   NWORDS  = THE TOTAL NUMBER OF WORDS FOUND IN CTEXT
C
C          (I*4)   LENTXT  = LENGTH IN BYTES OF 'CTEXT' STRING
C          (I*4)   IDELIM  = 0 => CTEXT CHARACTER IS NOT A DELIMITER
C                          > 0 => CTEXT CHARACTER IS A DELIMITER
C          (I*4)   ITOTAL  = NUMBER OF WORDS FOUND SO FAR
C          (I*4)   IINDEX  = IFIRST()/ILAST() INDEX OF CURRENT WORD
C          (I*4)   NLAST   = THE INDEX NO. OF THE LAST WORD TO EXTRACT
C          (I*4)   I       = GENERAL USE INDEX
C
C          (L*4)   LWORD   = .TRUE.  - PROCESSING AN IDENTIFIED WORD
C                            .FALSE. - PROCESSING SPACE BETWEEN WORDS
C
C ROUTINES: NONE
C
C NOTES:    IF THERE IS NO Nfirst WORD OR NO WORDS ARE FOUND
C           (I.E. INPUT STRING IS BLANK) THEN IWORDS=0
C
C AUTHOR:  PAUL E. BRIDEN (TESSELLA SUPPORT SERVICES PLC)
C          K1/0/37
C          JET EXT. 5023
C
C DATE:    20/05/93
C
C-----------------------------------------------------------------------
      INTEGER    NFIRST          , IWORDS           , NWORDS      ,
     &           LENTXT          , IDELIM           , ITOTAL      ,
     &           IINDEX          , NLAST            , I
C-----------------------------------------------------------------------
      LOGICAL    LWORD
C-----------------------------------------------------------------------
      CHARACTER  CTEXT*(*)       , CDELIM*(*)
C-----------------------------------------------------------------------
      INTEGER    IFIRST(IWORDS)  , ILAST(IWORDS)
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C
      LENTXT = LEN(CTEXT)
      NLAST  = IWORDS + NFIRST - 1
C
C-----------------------------------------------------------------------
C FIND THE REQUIRED WORDS
C-----------------------------------------------------------------------
C
         DO 1 I = 1,IWORDS
            IFIRST(I) = 0
            ILAST(I)  = 0
    1    CONTINUE
C
      ITOTAL = 0
      LWORD  = .FALSE.
      IINDEX = 0
C
         DO 2 I = 1,LENTXT

            IDELIM = INDEX( CDELIM , CTEXT(I:I) )

               IF (LWORD) THEN

                  IF (IDELIM.GT.0) THEN
                     LWORD  = .FALSE.
                     IF ((ITOTAL.GE.NFIRST).AND.(ITOTAL.LE.NLAST)) THEN
                        ILAST(IINDEX)  = I - 1
                     ENDIF
                  ENDIF

               ELSE

                  IF (IDELIM.EQ.0) THEN
                     LWORD  = .TRUE.
                     ITOTAL = ITOTAL + 1
                     IF ((ITOTAL.GE.NFIRST).AND.(ITOTAL.LE.NLAST)) THEN
                        IINDEX         = IINDEX + 1
                        IFIRST(IINDEX) = I
                        ILAST(IINDEX)  = LENTXT
                     ENDIF
                  ENDIF

               ENDIF

    2    CONTINUE
C
      IWORDS = IINDEX
      NWORDS = ITOTAL
C
C-----------------------------------------------------------------------
C
      RETURN
      END SUBROUTINE XXWORD
