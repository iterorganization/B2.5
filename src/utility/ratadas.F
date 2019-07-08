c     -------------------------------------------------------------------------
c     subroutine ratadas(nput,nc,rtt,rtnt,rtn,rtnn,
c                        rtsa,rtra,rtqa,rtcx,rtrd,rtbr,rtzz,rtz2,rtpt,rtpi,
c                        iyear,fix_recomb,ncnct,izpmin,izpmax)
c     -------------------------------------------------------------------------
c     INPUT: integer nput:     unit number for b2ar.prt output
c            integer nc:       nuclear charge
c            real (kind=R8) :: rtt(0:rtnt): temperatures in ev
c            integer rtnt:     number of temperatures
c            real (kind=R8) :: rtn(0:rtnn): densities in m**-3
c            integer rtnn:     number of densities
c            integer fix_recomb apply fix to include PRB
c     OUTPUT: real*8 rtsa(0:rtnt,0:rtnn,0:nc): ionisation r.c. in m**3/s
c             real (kind=R8) :: rtra(0:rtnt,0:rtnn,0:nc): recombination r.c. in m**3/s
c             real (kind=R8) :: rtqa(0:rtnt,0:rtnn,0:nc): electron heat loss r.c.
c                                            in eV*m**3/s
c             real (kind=R8) :: rtcx(0:rtnt,0:rtnn,0:nc): charge exchange r.c. in m**3/s
c             real (kind=R8) :: rtrd(0:rtnt,0:rtnn,0:nc): line radiation r.c.
c                                            in eV*m**3/s
c             real (kind=R8) :: rtbr(0:rtnt,0:rtnn,0:nc): bremsstrahlung radiation r.c.
c                                            in eV*m**3/s
c             real (kind=R8) :: rtzz(0:rtnt,0:rtnn,0:nc): effective charge in e
c             real (kind=R8) :: rtz2(0:rtnt,0:rtnn,0:nc): effective square charge in e**2
c             real (kind=R8) :: rtpt(0:rtnt,0:rtnn,0:nc): cumulative ionisation potential
c                                            in eV
c             real (kind=R8) :: rtpi(0:rtnt,0:rtnn,0:nc): effective ionisation potential
c                                            in eV
c             integer iyear     ADAS year of files
c                               (-1 if files comes from different years)
c             integer ncnct     number of partitions
c             integer izpmin()  array of lowest charge state of partition
c             integer izpmax()  array of highest charge state of partition
c     -------------------------------------------------------------------------
c     The paths to the ADAS files are taken from the file ratadas.filelist 
c     -------------------------------------------------------------------------
c     970826 bhs
c     -------------------------------------------------------------------------

      subroutine ratadas(nput,nc,rtt,rtnt,rtn,rtnn,
     &                   rtsa,rtra,rtqa,rtcx,rtrd,rtbr,
     &                   rtzz,rtz2,rtpt,rtpi,
     &                   iyear,fix_recomb,ncnct,izpmin,izpmax)
      use b2mod_types
      use b2mod_constants
      use b2mod_solpstop

      implicit none

c     INPUT
      integer nput
      integer nc
      integer rtnt
      integer rtnn
      integer fix_recomb
      real (kind=R8) ::
     *   rtt(0:rtnt)
      real (kind=R8) ::
     *   rtn(0:rtnn)
c     OUTPUT:
      real (kind=R8) ::
     *   rtsa(0:rtnt,0:rtnn,0:nc)
      real (kind=R8) ::
     *   rtra(0:rtnt,0:rtnn,0:nc)
      real (kind=R8) ::
     *   rtqa(0:rtnt,0:rtnn,0:nc)
      real (kind=R8) ::
     *   rtcx(0:rtnt,0:rtnn,0:nc)
      real (kind=R8) ::
     *   rtrd(0:rtnt,0:rtnn,0:nc)
      real (kind=R8) ::
     *   rtbr(0:rtnt,0:rtnn,0:nc)
      real (kind=R8) ::
     *   rtzz(0:rtnt,0:rtnn,0:nc)
      real (kind=R8) ::
     *   rtz2(0:rtnt,0:rtnn,0:nc)
      real (kind=R8) ::
     *   rtpt(0:rtnt,0:rtnn,0:nc)
      real (kind=R8) ::
     *   rtpi(0:rtnt,0:rtnn,0:nc)
      integer iyear, ncnct
      integer izpmin(0:nc), izpmax(0:nc)
c     LOCAL VARIABLES:
c     counter
      integer cc
c     filename acd (recombination)
      character fnrec*200
c     filename scd (ionisation)
      character fnion*200
c     filename ccd (charge exchange)
      character fncex*200
c     filename prb (bremsstrahlung + radiative recombination contribution)
      character fnprb*200
c     filename plt (line radiation)
      character fnplt*200
c     filename zcd (effective charge)
      character fnzav*200
c     filename ycd (effective square charge)
      character fnza2*200
c     filename ecd (effective ionisation potential)
      character fnpot*200
c     filename of actual file opened
      character filename*256, actual_file*256
c     helpers for composing filenames
      integer ibundling, iquote, quote_count
      character UIDIN*256, UIDSTD*256, ZEILE*80
      character SUBDIRrec*4,SUBDIRion*4,SUBDIRcex*4
      character SUBDIRprb*4,SUBDIRplt*4,SUBDIRzav*4
      character SUBDIRza2*4,SUBDIRpot*4
      character NEWDIRrec*4,NEWDIRion*4,NEWDIRcex*4
      character NEWDIRprb*4,NEWDIRplt*4,NEWDIRzav*4
      character NEWDIRza2*4,NEWDIRpot*4
      character ES*2, EB*2, Bdescr*8
c     helpers for reading in data
      integer IZMAX,IZ1MIN,IZ1MAX
      integer isdimd,iddimd,itdimd
      parameter (isdimd=92+1,iddimd=100,itdimd=100)
      integer nptnl,ndptnl,ndptn,ndptnc,ndcnct
      parameter (ndptnl=4,ndptn=128,ndptnc=256,ndcnct=100)
      integer nptn(ndptnl),nptnc(ndptnl,ndptn)
      integer iptnla(ndptnl),iptna(ndptnl,ndptn)
      integer iptnca(ndptnl,ndptn,ndptnc)
      integer icnctv(ndcnct)
      integer iblmx,ismax
      integer isppr(isdimd),ispbr(isdimd),isstgr(isdimd)
      integer IDMAXDrec,ITMAXDrec,IDMAXDion,ITMAXDion,
     &        IDMAXDcex,ITMAXDcex,
     &        IDMAXDprb,ITMAXDprb,
     &        IDMAXDplt,ITMAXDplt,
     &        IDMAXDzav,ITMAXDzav,
     &        IDMAXDza2,ITMAXDza2,
     &        IDMAXDpot,ITMAXDpot
      logical lres,lstan,lptn,ex
      real (kind=R8) :: dnr_ams
      character dnr_ele*12
c     read in data
      integer maxne,ic
c     log10 of densities (in 1/cm3), temperatures (in eV)
      real (kind=R8) ::
     *   DDENSDrec(iddimd),DTEVDrec(itdimd)
      real (kind=R8) ::
     *   DDENSDion(iddimd),DTEVDion(itdimd)
      real (kind=R8) ::
     *   DDENSDcex(iddimd),DTEVDcex(itdimd)
      real (kind=R8) ::
     *   DDENSDprb(iddimd),DTEVDprb(itdimd)
      real (kind=R8) ::
     *   DDENSDplt(iddimd),DTEVDplt(itdimd)
      real (kind=R8) ::
     *   DDENSDzav(iddimd),DTEVDzav(itdimd)
      real (kind=R8) ::
     *   DDENSDza2(iddimd),DTEVDza2(itdimd)
      real (kind=R8) ::
     *   DDENSDpot(iddimd),DTEVDpot(itdimd)
      real (kind=R8) ::
     *   Ddummy(iddimd),Tdummy(itdimd)
c     log10 of rate coefficients (in cm3/s)
      real (kind=R8) ::
     *   DREC(isdimd,itdimd,iddimd),DION(isdimd,itdimd,iddimd),
     &   DCEX(isdimd,itdimd,iddimd),DPLT(isdimd,itdimd,iddimd),
     &   DPRB(isdimd,itdimd,iddimd),DTMP1(itdimd,iddimd),
     1   DTMP2(itdimd,iddimd),DTMP3(itdimd,iddimd)
c     effective coefficients
      real (kind=R8) ::
     *   DZAV(isdimd,itdimd,iddimd),DZA2(isdimd,itdimd,iddimd),
     &   DPOT(isdimd,itdimd,iddimd)
c     helpers for filling ZTER,ZNER,rtsa,rtra,rtqq,rtcx
c     counter for temperature,density
      integer cT,cD
c     log10 of temperature,density
      real (kind=R8) ::
     *   logT,logD
c     return value
      real (kind=R8) ::
     *   ans
c     ionisation potentials
      real (kind=R8) ::
     *   ip(0:nc-1)

      logical silent_open
      external silent_open, subini, subend, xerrab

*   ..subprogram start-up calls
      call subini ('ratadas')

      rtsa=0.0_R8
      rtra=0.0_R8
      rtqa=0.0_R8
      rtcx=0.0_R8
      rtrd=0.0_R8
      rtbr=0.0_R8
      rtzz=-1.0_R8
      rtz2=-1.0_R8
      rtpi=-1.0_R8
      if(fix_recomb.ne.0) then 
        write(*,*) 'RATADAS: FIX_RECOMB <> 0 : ',fix_recomb 
      endif 
      do ic = 1, itdimd
        Tdummy(ic) = (-1.0_R8+5.0_R8*(ic-1.0_R8)/(itdimd-1.0_R8))
      enddo
      do ic = 1, iddimd
        Ddummy(ic) = (12.0_R8+9.0_R8*(ic-1.0_R8)/(iddimd-1.0_R8))
      enddo

c     reading in filename and mass
      filename='ratadas.filelist'
      uidstd=trim(solpstop()) // '/modules/adas'
      call find_file(filename,ex)
      call xertst(silent_open(99,filename,
     1 'data.local|modules/B2.5/Database',actual_file), 'missing file')
      write(*,'(a/a)')
     1 'Getting ADAS file list from ',trim(actual_file)
      write(nput,'(a/a)')
     1 'Getting ADAS file list from ',trim(actual_file)
      read(99,*) UIDIN
      if(uidin(1:8).eq.'SOLPSTOP') then
        if (uidin(9:9).ne.'/') then
          uidin=trim(solpstop()) // '/' // trim(uidin(9:))
        else
          uidin=trim(solpstop()) // '/' // trim(uidin(10:))
        endif
        write(*,'(a,a)') 'ADAS UIDIN changed to ',trim(UIDIN)
      endif
      read(99,*) maxne
      if (nc.gt.maxne) then
         write(*,*) 'ADAStables: element number',nc,'not supported!'
         return
      endif
c     jump over the preceding lines
      do cc=1,nc-1
         read(99,*) SUBDIRrec,SUBDIRion,SUBDIRcex,SUBDIRprb,
     &              SUBDIRplt,SUBDIRzav,SUBDIRpot,ES
      enddo
c     now read the relevant line
      read(99,'(a80)') ZEILE
      quote_count = 0
      iquote = index(ZEILE,'''')
      do while (iquote.gt.0)
         quote_count=quote_count+1
         if (index(ZEILE(iquote+1:80),'''').ne.0) then
            iquote = index(ZEILE(iquote+1:80),'''') + iquote
         else
            iquote = 0
         endif
      enddo
      if (quote_count.eq.16) then
         read(ZEILE,*) SUBDIRrec,SUBDIRion,SUBDIRcex,SUBDIRprb,
     &                 SUBDIRplt,SUBDIRzav,SUBDIRpot,ES
      else if (quote_count.eq.18) then
         read(ZEILE,*) SUBDIRrec,SUBDIRion,SUBDIRcex,SUBDIRprb,
     &                 SUBDIRplt,SUBDIRzav,SUBDIRpot,ES,Bdescr
      else
         call xerrab ('Incomplete line found in ratadas.filelist!')
      endif
      SUBDIRza2=SUBDIRzav
      close(99)

c     check ADAS year
      if((SUBDIRrec.ne.'none') .and.
     &   (SUBDIRrec.eq.SUBDIRion .or. SUBDIRion.eq.'none') .and.
     &   (SUBDIRrec.eq.SUBDIRcex .or. SUBDIRcex.eq.'none') .and.
     &   (SUBDIRrec.eq.SUBDIRprb .or. SUBDIRprb.eq.'none') .and.
     &   (SUBDIRrec.eq.SUBDIRplt .or. SUBDIRplt.eq.'none') .and.
     &   (SUBDIRrec.eq.SUBDIRzav .or. SUBDIRzav.eq.'none') .and.
     &   (SUBDIRrec.eq.SUBDIRpot .or. SUBDIRpot.eq.'none') ) then
          read(SUBDIRrec,*) iyear
          if (iyear.gt.100) then
            ibundling = mod(iyear,100)
            iyear = iyear/100
            if (ibundling.le.9) then
              write(EB,'(a1,i1)') '0',ibundling
            else
              write(EB,'(i2)') ibundling
            endif
          else
            ibundling = 0
          endif
          write(NEWDIRrec,'(i2,2x)') iyear
          write(NEWDIRion,'(i2,2x)') iyear
          write(NEWDIRcex,'(i2,2x)') iyear
          write(NEWDIRprb,'(i2,2x)') iyear
          write(NEWDIRplt,'(i2,2x)') iyear
          write(NEWDIRzav,'(i2,2x)') iyear
          write(NEWDIRza2,'(i2,2x)') iyear
          write(NEWDIRpot,'(i2,2x)') iyear
      else
          ibundling = 0
          iyear = -1
          NEWDIRrec = SUBDIRrec
          NEWDIRion = SUBDIRion
          NEWDIRcex = SUBDIRcex
          NEWDIRprb = SUBDIRprb
          NEWDIRplt = SUBDIRplt
          NEWDIRzav = SUBDIRzav
          NEWDIRza2 = SUBDIRza2
          NEWDIRpot = SUBDIRpot
      endif

c     ADAS file names
      if (ibundling.eq.0) then
        fnrec='/adf11/acd'//NEWDIRrec//'/acd'//NEWDIRrec//
     &      '_'//ES//'.dat'
        fnion='/adf11/scd'//NEWDIRion//'/scd'//NEWDIRion//
     &      '_'//ES//'.dat'
        fncex='/adf11/ccd'//NEWDIRcex//'/ccd'//NEWDIRcex//
     &      '_'//ES//'.dat'
        fnprb='/adf11/prb'//NEWDIRprb//'/prb'//NEWDIRprb//
     &      '_'//ES//'.dat'
        fnplt='/adf11/plt'//NEWDIRplt//'/plt'//NEWDIRplt//
     &      '_'//ES//'.dat'
        fnzav='/adf11/zcd'//NEWDIRzav//'/zcd'//NEWDIRzav//
     &      '_'//ES//'.dat'
        fnza2='/adf11/ycd'//NEWDIRza2//'/ycd'//NEWDIRza2//
     &      '_'//ES//'.dat'
        fnpot='/adf11/ecd'//NEWDIRpot//'/ecd'//NEWDIRpot//
     &      '_'//ES//'.dat'
      else if (quote_count.eq.16) then
        fnrec='/adf11/acd'//NEWDIRrec//'/acd'//NEWDIRrec//
     &      '_'//ES//'_'//EB//'.dat'
        fnion='/adf11/scd'//NEWDIRion//'/scd'//NEWDIRion//
     &      '_'//ES//'_'//EB//'.dat'
        fncex='/adf11/ccd'//NEWDIRcex//'/ccd'//NEWDIRcex//
     &      '_'//ES//'_'//EB//'.dat'
        fnprb='/adf11/prb'//NEWDIRprb//'/prb'//NEWDIRprb//
     &      '_'//ES//'_'//EB//'.dat'
        fnplt='/adf11/plt'//NEWDIRplt//'/plt'//NEWDIRplt//
     &      '_'//ES//'_'//EB//'.dat'
        fnzav='/adf11/zcd'//NEWDIRzav//'/zcd'//NEWDIRzav//
     &      '_'//ES//'_'//EB//'.dat'
        fnza2='/adf11/ycd'//NEWDIRza2//'/ycd'//NEWDIRza2//
     &      '_'//ES//'_'//EB//'.dat'
        fnpot='/adf11/ecd'//NEWDIRpot//'/ecd'//NEWDIRpot//
     &      '_'//ES//'_'//EB//'.dat'
      else if (quote_count.eq.18) then
        fnrec='/adf11/acd'//NEWDIRrec//'/acd'//NEWDIRrec//
     &      '_'//ES//'_'//EB//'_'//Bdescr//'.dat'
        fnion='/adf11/scd'//NEWDIRion//'/scd'//NEWDIRion//
     &      '_'//ES//'_'//EB//'_'//Bdescr//'.dat'
        fncex='/adf11/ccd'//NEWDIRcex//'/ccd'//NEWDIRcex//
     &      '_'//ES//'_'//EB//'_'//Bdescr//'.dat'
        fnprb='/adf11/prb'//NEWDIRprb//'/prb'//NEWDIRprb//
     &      '_'//ES//'_'//EB//'_'//Bdescr//'.dat'
        fnplt='/adf11/plt'//NEWDIRplt//'/plt'//NEWDIRplt//
     &      '_'//ES//'_'//EB//'_'//Bdescr//'.dat'
        fnzav='/adf11/zcd'//NEWDIRzav//'/zcd'//NEWDIRzav//
     &      '_'//ES//'_'//EB//'_'//Bdescr//'.dat'
        fnza2='/adf11/ycd'//NEWDIRza2//'/ycd'//NEWDIRza2//
     &      '_'//ES//'_'//EB//'_'//Bdescr//'.dat'
        fnpot='/adf11/ecd'//NEWDIRpot//'/ecd'//NEWDIRpot//
     &      '_'//ES//'_'//EB//'_'//Bdescr//'.dat'
      endif
      call strip_spaces(fnrec)
      call strip_spaces(fnion)
      call strip_spaces(fncex)
      call strip_spaces(fnprb)
      call strip_spaces(fnplt)
      call strip_spaces(fnzav)
      call strip_spaces(fnza2)
      call strip_spaces(fnpot)
      filename=trim(UIDIN)//fnrec
      call find_file(filename,ex)
      if(ex) then
        fnrec=trim(filename)
      else
        fnrec=trim(UIDSTD)//fnrec
      endif
      filename=trim(UIDIN)//fnion
      call find_file(filename,ex)
      if(ex) then
        fnion=trim(filename)
      else
        fnion=trim(UIDSTD)//fnion
      endif
      filename=trim(UIDIN)//fncex
      call find_file(filename,ex)
      if(ex) then
        fncex=trim(filename)
      else
        fncex=trim(UIDSTD)//fncex
      endif
      filename=trim(UIDIN)//fnprb
      call find_file(filename,ex)
      if(ex) then
        fnprb=trim(filename)
      else
        fnprb=trim(UIDSTD)//fnprb
      endif
      filename=trim(UIDIN)//fnplt
      call find_file(filename,ex)
      if(ex) then
        fnplt=trim(filename)
      else
        fnplt=trim(UIDSTD)//fnplt
      endif
      filename=trim(UIDIN)//fnzav
      call find_file(filename,ex)
      if(ex) then
        fnzav=trim(filename)
      else
        fnzav=trim(UIDSTD)//fnzav
      endif
      filename=trim(UIDIN)//fnza2
      call find_file(filename,ex)
      if(ex) then
        fnza2=trim(filename)
      else
        fnza2=trim(UIDSTD)//fnza2
      endif
      filename=trim(UIDIN)//fnpot
      call find_file(filename,ex)
      if(ex) then
        fnpot=trim(filename)
      else
        fnpot=trim(UIDSTD)//fnpot
      endif

      write(*,*) 'Using ADAS files:'
      if (SUBDIRrec.ne.'none') then
        write(*,'(a)') trim(fnrec)
      else
        write(*,*) 'Skipping recombination data'
      endif
      if (SUBDIRion.ne.'none') then
        write(*,'(a)') trim(fnion)
      else
        write(*,*) 'Skipping ionisation data'
      endif
      if (SUBDIRcex.ne.'none') then
        write(*,'(a)') trim(fncex)
      else
        write(*,*) 'Skipping charge exchange data'
      endif
      if (SUBDIRprb.ne.'none') then
        write(*,'(a)') trim(fnprb)
      else
        write(*,*) 'Skipping continuum radiation data'
      endif
      if (SUBDIRplt.ne.'none') then
        write(*,'(a)') trim(fnplt)
      else
        write(*,*) 'Skipping line radiation data'
      endif
      if (SUBDIRzav.ne.'none') then
        write(*,'(a)') trim(fnzav)
      else
        write(*,*) 'Skipping effective charge data'
      endif
      if (SUBDIRza2.ne.'none') then
        write(*,'(a)') trim(fnza2)
      else
        write(*,*) 'Skipping effective square charge data'
      endif
      if (SUBDIRpot.ne.'none') then
        write(*,'(a)') trim(fnpot)
      else
        write(*,*) 'Skipping effective ionisation potential data'
      endif


c     reading in data, recombination
      if(SUBDIRrec.eq.'none') then
        write(*,*) 
     .   'Setting recombination to zero'
        write(nput,*) 
     .   'Setting recombination to zero'
        IDMAXDrec=iddimd
        ITMAXDrec=itdimd
        DDENSDrec(1:IDMAXDrec)=Ddummy(1:IDMAXDrec)
        DTEVDrec(1:ITMAXDrec)=Tdummy(1:ITMAXDrec)
        DREC=-74.0_R8
      else 
     . if(silent_open(99,fnrec,'data.local|modules',actual_file)) then
        write(*,'(a,a,a)') 'Used ',trim(actual_file),
     .   ' for recombination'
        write(nput,'(a,a,a)') 'Used ',trim(actual_file),
     .   ' for recombination'
        call xxdata_11(99,1,
     *           isdimd,iddimd,itdimd,ndptnl,ndptn,ndptnc,ndcnct,
     *           IZMAX,IZ1MIN,IZ1MAX,
     *           nptnl,nptn,nptnc,iptnla,iptna,iptnca,
     *           ncnct,icnctv,iblmx,ismax,dnr_ele,dnr_ams,
     *           isppr,ispbr,isstgr,
     *           IDMAXDrec,ITMAXDrec,DDENSDrec,DTEVDrec,DREC,
     *           lres,lstan,lptn)
        if (IZMAX.ne.nc) then
          write(*,*) 'ADAS tables: wrong element.'
          write(*,*) 'reading file: ',fnrec
          return
        endif
        if ((IDMAXDrec.gt.iddimd).or.(ITMAXDrec.gt.itdimd)) then
          write(*,*) 'ADAS tables: dimensions too big.'
          write(*,*) 'Change iddimd and itdimd in source code.'
          write(*,*) 'reading file: ',fnrec
          return
        endif
        close(99)
      else
        goto 990
      endif
c     reading in data, ionisation
      if(SUBDIRion.eq.'none') then
        write(*,*) 
     .   'Setting ionisation to zero'
        write(nput,*) 
     .   'Setting ionisation to zero'
        IDMAXDion=iddimd
        ITMAXDion=itdimd
        DDENSDion(1:IDMAXDion)=Ddummy(1:IDMAXDion)
        DTEVDion(1:ITMAXDion)=Tdummy(1:ITMAXDion)
        DION=-74.0_R8
      else 
     . if(silent_open(99,fnion,'data.local|modules',actual_file)) then
        write(*,'(a,a,a)') 'Used ',trim(actual_file),
     .   ' for ionisation'
        write(nput,'(a,a,a)') 'Used ',trim(actual_file),
     .   ' for ionisation'
        call xxdata_11(99,2,
     *           isdimd,iddimd,itdimd,ndptnl,ndptn,ndptnc,ndcnct,
     *           IZMAX,IZ1MIN,IZ1MAX,
     *           nptnl,nptn,nptnc,iptnla,iptna,iptnca,
     *           ncnct,icnctv,iblmx,ismax,dnr_ele,dnr_ams,
     *           isppr,ispbr,isstgr,
     *           IDMAXDion,ITMAXDion,DDENSDion,DTEVDion,DION,
     *           lres,lstan,lptn)
        if (IZMAX.ne.nc) then
          write(*,*) 'ADAS tables: wrong element.'
          write(*,*) 'reading file: ',fnion
          return
        endif
        if ((IDMAXDion.gt.iddimd).or.(ITMAXDion.gt.itdimd)) then
          write(*,*) 'ADAS tables: dimensions too big.'
          write(*,*) 'Change iddimd and itdimd in source code.'
          write(*,*) 'reading file: ',fnion
          return
        endif
        close(99) 
      else
        goto 990
      endif
c     reading in data, charge exchange
      if(SUBDIRcex.eq.'none') then
        write(*,*) 
     .   'Setting charge exchange to zero'
        write(nput,*) 
     .   'Setting charge exchange to zero'
        IDMAXDcex=iddimd
        ITMAXDcex=itdimd
        DDENSDcex(1:IDMAXDcex)=Ddummy(1:IDMAXDcex)
        DTEVDcex(1:ITMAXDcex)=Tdummy(1:ITMAXDcex)
        DCEX=-74.0_R8
      else 
     . if(silent_open(99,fncex,'data.local|modules',actual_file)) then
        write(*,'(a,a,a)') 'Used ',trim(actual_file),
     .   ' for charge exchange'
        write(nput,'(a,a,a)') 'Used ',trim(actual_file),
     .   ' for charge exchange'
        call xxdata_11(99,3,
     *           isdimd,iddimd,itdimd,ndptnl,ndptn,ndptnc,ndcnct,
     *           IZMAX,IZ1MIN,IZ1MAX,
     *           nptnl,nptn,nptnc,iptnla,iptna,iptnca,
     *           ncnct,icnctv,iblmx,ismax,dnr_ele,dnr_ams,
     *           isppr,ispbr,isstgr,
     *           IDMAXDcex,ITMAXDcex,DDENSDcex,DTEVDcex,DCEX,
     *           lres,lstan,lptn)
        if (IZMAX.ne.nc) then
          write(*,*) 'ADAS tables: wrong element.'
          write(*,*) 'reading file: ',fncex
          return
        endif
        if ((IDMAXDcex.gt.iddimd).or.(ITMAXDcex.gt.itdimd)) then
          write(*,*) 'ADAS tables: dimensions too big.'
          write(*,*) 'Change iddimd and itdimd in source code.'
          write(*,*) 'reading file: ',fncex
          return
        endif
        close(99) 
      else
        goto 990
      endif
c     reading in data, bremsstrahlung + radiative recombination contribution
      if(SUBDIRprb.eq.'none') then
        write(*,*) 
     .   'Setting continuum radiation to zero'
        write(nput,*) 
     .   'Setting continuum radiation to zero'
        IDMAXDprb=iddimd
        ITMAXDprb=itdimd
        DDENSDprb(1:IDMAXDprb)=Ddummy(1:IDMAXDprb)
        DTEVDprb(1:ITMAXDprb)=Tdummy(1:ITMAXDprb)
        DPRB=-74.0_R8
      else 
     . if(silent_open(99,fnprb,'data.local|modules',actual_file)) then
        write(*,'(a,a,a)') 'Used ',trim(actual_file),
     .   ' for bremsstrahlung'
        write(nput,'(a,a,a)') 'Used ',trim(actual_file),
     .   ' for bremsstrahlung'
        call xxdata_11(99,4,
     *           isdimd,iddimd,itdimd,ndptnl,ndptn,ndptnc,ndcnct,
     *           IZMAX,IZ1MIN,IZ1MAX,
     *           nptnl,nptn,nptnc,iptnla,iptna,iptnca,
     *           ncnct,icnctv,iblmx,ismax,dnr_ele,dnr_ams,
     *           isppr,ispbr,isstgr,
     *           IDMAXDprb,ITMAXDprb,DDENSDprb,DTEVDprb,DPRB,
     *           lres,lstan,lptn)
        if (IZMAX.ne.nc) then
          write(*,*) 'ADAS tables: wrong element.'
          write(*,*) 'reading file: ',fnprb
          return
        endif
        if ((IDMAXDprb.gt.iddimd).or.(ITMAXDprb.gt.itdimd)) then
          write(*,*) 'ADAS tables: dimensions too big.'
          write(*,*) 'Change iddimd and itdimd in source code.'
          write(*,*) 'reading file: ',fnprb
          return
        endif
        close(99) 
      else
        goto 990
      endif
c     reading in data, line radiation
      if(SUBDIRplt.eq.'none') then
        write(*,*) 
     .   'Setting line radiation to zero'
        write(nput,*) 
     .   'Setting line radiation to zero'
        IDMAXDplt=iddimd
        ITMAXDplt=itdimd
        DDENSDplt(1:IDMAXDplt)=Ddummy(1:IDMAXDplt)
        DTEVDplt(1:ITMAXDplt)=Tdummy(1:ITMAXDplt)
        DPLT=-74.0_R8
      else 
     . if(silent_open(99,fnplt,'data.local|modules',actual_file)) then
        write(*,'(a,a,a)') 'Used ',trim(actual_file),
     .   ' for line radiation'
        write(nput,'(a,a,a)') 'Used ',trim(actual_file),
     .   ' for line radiation'
        call xxdata_11(99,8,
     *           isdimd,iddimd,itdimd,ndptnl,ndptn,ndptnc,ndcnct,
     *           IZMAX,IZ1MIN,IZ1MAX,
     *           nptnl,nptn,nptnc,iptnla,iptna,iptnca,
     *           ncnct,icnctv,iblmx,ismax,dnr_ele,dnr_ams,
     *           isppr,ispbr,isstgr,
     *           IDMAXDplt,ITMAXDplt,DDENSDplt,DTEVDplt,DPLT,
     *           lres,lstan,lptn)
        if (IZMAX.ne.nc) then
          write(*,*) 'ADAS tables: wrong element.'
          write(*,*) 'reading file: ',fnplt
          return
        endif
        if ((IDMAXDplt.gt.iddimd).or.(ITMAXDplt.gt.itdimd)) then
          write(*,*) 'ADAS tables: dimensions too big.'
          write(*,*) 'Change iddimd and itdimd in source code.'
          write(*,*) 'reading file: ',fnplt
          return
        endif
        close(99) 
      else
        goto 990
      endif
c     reading in data, effective charge
      if(SUBDIRzav.eq.'none') then
        write(*,*) 
     .   'Setting default effective charge'
        write(nput,*) 
     .   'Setting default effective charge'
        IDMAXDzav=iddimd
        ITMAXDzav=itdimd
        DDENSDzav(1:IDMAXDzav)=Ddummy(1:IDMAXDzav)
        DTEVDzav(1:ITMAXDzav)=Tdummy(1:ITMAXDzav)
        DZAV=-1.0_R8
        ncnct=nc+1
        do ic = 0, nc
          izpmin(ic) = ic
          izpmax(ic) = ic
        enddo
      else 
     . if(silent_open(99,fnzav,'data.local|modules',actual_file)) then
        write(*,'(a,a,a)') 'Used ',trim(actual_file),
     .   ' for effective charge'
        write(nput,'(a,a,a)') 'Used ',trim(actual_file),
     .   ' for effective charge'
        call xxdata_11(99,10,
     *           isdimd,iddimd,itdimd,ndptnl,ndptn,ndptnc,ndcnct,
     *           IZMAX,IZ1MIN,IZ1MAX,
     *           nptnl,nptn,nptnc,iptnla,iptna,iptnca,
     *           ncnct,icnctv,iblmx,ismax,dnr_ele,dnr_ams,
     *           isppr,ispbr,isstgr,
     *           IDMAXDzav,ITMAXDzav,DDENSDzav,DTEVDzav,DZAV,
     *           lres,lstan,lptn)
        if (IZMAX.ne.nc) then
          write(*,*) 'ADAS tables: wrong element.'
          write(*,*) 'reading file: ',fnzav
          return
        endif
        if ((IDMAXDzav.gt.iddimd).or.(ITMAXDzav.gt.itdimd)) then
          write(*,*) 'ADAS tables: dimensions too big.'
          write(*,*) 'Change iddimd and itdimd in source code.'
          write(*,*) 'reading file: ',fnzav
          return
        endif
        close(99)
        do ic = 1, ncnct
          izpmin(ic-1) = iptnca(1,ic,1)
          izpmax(ic-1) = iptnca(1,ic,nptnc(1,ic))
        enddo 
      else
        write(*,*) 
     .   'Using default values for effective charge'
        write(nput,*) 
     .   'Using default values for effective charge'
        IDMAXDzav=IDMAXDplt
        ITMAXDzav=ITMAXDplt
        DDENSDzav(1:IDMAXDplt)=DDENSDplt(1:IDMAXDplt)
        DTEVDzav(1:ITMAXDplt)=DTEVDplt(1:ITMAXDplt)
        DZAV=-1.0_R8
        ncnct=nc+1
        do ic = 0, nc
          izpmin(ic) = ic
          izpmax(ic) = ic
        enddo
      endif
c     reading in data, effective square charge
      if(SUBDIRza2.eq.'none') then
        write(*,*) 
     .   'Setting default effective square charge'
        write(nput,*) 
     .   'Setting default effective square charge'
        IDMAXDza2=iddimd
        ITMAXDza2=itdimd
        DDENSDza2(1:IDMAXDza2)=Ddummy(1:IDMAXDza2)
        DTEVDza2(1:ITMAXDza2)=Tdummy(1:ITMAXDza2)
        DZA2=-1.0_R8
      else 
     . if(silent_open(99,fnza2,'data.local|modules',actual_file)) then
        write(*,'(a,a,a)') 'Used ',trim(actual_file),
     .   ' for effective square charge'
        write(nput,'(a,a,a)') 'Used ',trim(actual_file),
     .   ' for effective square charge'
        call xxdata_11(99,11,
     *           isdimd,iddimd,itdimd,ndptnl,ndptn,ndptnc,ndcnct,
     *           IZMAX,IZ1MIN,IZ1MAX,
     *           nptnl,nptn,nptnc,iptnla,iptna,iptnca,
     *           ncnct,icnctv,iblmx,ismax,dnr_ele,dnr_ams,
     *           isppr,ispbr,isstgr,
     *           IDMAXDza2,ITMAXDza2,DDENSDza2,DTEVDza2,DZA2,
     *           lres,lstan,lptn)
        if (IZMAX.ne.nc) then
          write(*,*) 'ADAS tables: wrong element.'
          write(*,*) 'reading file: ',fnza2
          return
        endif
        if ((IDMAXDza2.gt.iddimd).or.(ITMAXDza2.gt.itdimd)) then
          write(*,*) 'ADAS tables: dimensions too big.'
          write(*,*) 'Change iddimd and itdimd in source code.'
          write(*,*) 'reading file: ',fnza2
          return
        endif
        close(99) 
      else
        write(*,*) 
     .   'Using default values for effective square charge'
        write(nput,*) 
     .   'Using default values for effective square charge'
        IDMAXDza2=IDMAXDplt
        ITMAXDza2=ITMAXDplt
        DDENSDza2(1:IDMAXDplt)=DDENSDplt(1:IDMAXDplt)
        DTEVDza2(1:ITMAXDplt)=DTEVDplt(1:ITMAXDplt)
        DZA2=-1.0_R8
      endif
c     reading in data, effective ionisation potential
      if(SUBDIRpot.eq.'none') then
        write(*,*) 
     .   'Setting default effective ionisation potential'
        write(nput,*) 
     .   'Setting default effective ionisation potential'
        IDMAXDpot=iddimd
        ITMAXDpot=itdimd
        DDENSDpot(1:IDMAXDpot)=Ddummy(1:IDMAXDpot)
        DTEVDpot(1:ITMAXDpot)=Tdummy(1:ITMAXDpot)
        DPOT=-1.0_R8
      else 
     . if(silent_open(99,fnpot,'data.local|modules',actual_file)) then
        write(*,'(a,a,a)') 'Used ',trim(actual_file),
     .   ' for effective ionisation potential'
        write(nput,'(a,a,a)') 'Used ',trim(actual_file),
     .   ' for effective ionisation potential'
        call xxdata_11(99,12,
     *           isdimd,iddimd,itdimd,ndptnl,ndptn,ndptnc,ndcnct,
     *           IZMAX,IZ1MIN,IZ1MAX,
     *           nptnl,nptn,nptnc,iptnla,iptna,iptnca,
     *           ncnct,icnctv,iblmx,ismax,dnr_ele,dnr_ams,
     *           isppr,ispbr,isstgr,
     *           IDMAXDpot,ITMAXDpot,DDENSDpot,DTEVDpot,DPOT,
     *           lres,lstan,lptn)
        if (IZMAX.ne.nc) then
          write(*,*) 'ADAS tables: wrong element.'
          write(*,*) 'reading file: ',fnpot
          return
        endif
        if ((IDMAXDpot.gt.iddimd).or.(ITMAXDpot.gt.itdimd)) then
          write(*,*) 'ADAS tables: dimensions too big.'
          write(*,*) 'Change iddimd and itdimd in source code.'
          write(*,*) 'reading file: ',fnpot
          return
        endif
        close(99) 
      else
        write(*,*) 
     .   'Using default values for effective ionisation potential'
        write(nput,*) 
     .   'Using default values for effective ionisation potential'
        IDMAXDpot=IDMAXDplt
        ITMAXDpot=ITMAXDplt
        DDENSDpot(1:IDMAXDplt)=DDENSDplt(1:IDMAXDplt)
        DTEVDpot(1:ITMAXDplt)=DTEVDplt(1:ITMAXDplt)
        DPOT=-1.0_R8
      endif
c     reading in data, ionisation potentials
      call xertst (silent_open(99,'ionization_potentials',
     1 'data.local|modules/B2.5/Database',actual_file), 'missing file')
      write(*,'(a,a,a)') 'Used ',trim(actual_file),
     .   ' for ionisation potentials table'
      write(nput,'(a,a,a)') 'Used ',trim(actual_file),
     .   ' for ionisation potentials table'
      do ic=1,nc-1
         read(99,*)
      enddo
      read(99,*) ip(0:nc-1)
      close(99)
c     filling rtsa(ionisation)
      write(*,*) 'rtsa'
      do ic=0,ncnct-2
        dtmp1(1:ITMAXDion,1:IDMAXDion)=
     .   DION(ic+1,1:ITMAXDion,1:IDMAXDion)
        do cT=0,rtnt
          logT=log10(rtt(cT))
          do cD=0,rtnn
            logD=log10(rtn(cD)/1.0e6_R8)
            call int2d(DTMP1(1:itdimd,1:iddimd),itdimd,iddimd,
     1       DTEVDion,ITMAXDion,DDENSDion,IDMAXDion,
     2       logT,logD,ans)
            rtsa(cT,cD,ic)=10.0_R8**(ans-6.0_R8)
          enddo
        enddo
      enddo
c     filling rtra(recombination)
      write(*,*) 'rtra'
      do ic=1,ncnct-1
        dtmp1(1:ITMAXDrec,1:IDMAXDrec)=
     .   DREC(ic,1:ITMAXDrec,1:IDMAXDrec)
        do cT=0,rtnt
          logT=log10(rtt(cT))
          do cD=0,rtnn
            logD=log10(rtn(cD)/1.0e6_R8)
            call int2d(DTMP1(1:itdimd,1:iddimd),itdimd,iddimd,
     1       DTEVDrec,ITMAXDrec,DDENSDrec,IDMAXDrec,
     2       logT,logD,ans)
            rtra(cT,cD,ic)=10.0_R8**(ans-6.0_R8)
          enddo
        enddo
      enddo
c     filling rtqa(electron heat loss), rtrd(line radiation) and rtbr(bremsstrahlung)
      write(*,*) 'rtqa,rtrd,rtbr'
      do ic=0,ncnct-1
        if(ic.gt.0) then
         dtmp1(1:ITMAXDprb,1:IDMAXDprb)=
     .    DPRB(ic,1:ITMAXDprb,1:IDMAXDprb)
        endif
         if(ic.lt.ncnct-1) then
          dtmp2(1:ITMAXDplt,1:IDMAXDplt)=
     .     DPLT(ic+1,1:ITMAXDplt,1:IDMAXDplt)
          dtmp3(1:ITMAXDpot,1:IDMAXDpot)=
     .     DPOT(ic+2,1:ITMAXDpot,1:IDMAXDpot)
        endif
        do cT=0,rtnt
          logT=log10(rtt(cT))
          do cD=0,rtnn
            logD=log10(rtn(cD)/1.0e6_R8)
            rtrd(cT,cD,ic)=0.0_R8
            rtbr(cT,cD,ic)=1.69e-38_R8/qe*ic**2*sqrt(rtt(cT))
            if(fix_recomb.eq.0) then 
              rtqa(cT,cD,ic)=rtbr(cT,cD,ic) 
            else 
              if(ic.gt.0) then 
                call int2d(DTMP1(1:itdimd,1:iddimd),itdimd,iddimd, 
     1           DTEVDprb,ITMAXDprb,DDENSDprb,IDMAXDprb, 
     2           logT,logD,ans) 
                rtqa(cT,cD,ic)=10.0_R8**(ans-6.0_R8)/qe 
                if(rtqa(cT,cD,ic).ge.rtbr(cT,cD,ic)) then 
                  rtrd(cT,cD,ic)=rtqa(cT,cD,ic)-rtbr(cT,cD,ic) 
                else 
                  rtbr(cT,cD,ic)=rtqa(cT,cD,ic) 
                  rtrd(cT,cD,ic)=0.0_R8
                endif 
c there should also be an electron heating term corresponding to 
c ip(ic-1)*rtra(cT,cD,ic) but this could make the cooling rate negative and 
c so is included explicitly in b2stel if fix_recomb is nonzero
              else 
                rtqa(cT,cD,ic)=0.0_R8 
              endif 
            endif 
            if(ic.lt.ncnct-1) then
              call int2d(DTMP2(1:itdimd,1:iddimd),itdimd,iddimd,
     1         DTEVDplt,ITMAXDplt,DDENSDplt,IDMAXDplt,
     2         logT,logD,ans)
              rtqa(cT,cD,ic)=rtqa(cT,cD,ic)+10.0_R8**(ans-6.0_R8)/qe
              rtrd(cT,cD,ic)=rtrd(cT,cD,ic)+10.0_R8**(ans-6.0_R8)/qe
              call int2d(DTMP3(1:itdimd,1:iddimd),itdimd,iddimd,
     1         DTEVDpot,ITMAXDpot,DDENSDpot,IDMAXDpot,
     2         logT,logD,ans)
              if (ans.lt.0.0_R8) then
                rtqa(cT,cD,ic)=rtqa(cT,cD,ic)+ip(ic)*rtsa(cT,cD,ic)
              else
                rtqa(cT,cD,ic)=rtqa(cT,cD,ic)+ans*rtsa(cT,cD,ic)
              endif
            endif
          enddo
        enddo
      enddo
c     filling rtcx(charge exchange)
      write(*,*) 'rtcx'
      do ic=1,ncnct-1
        dtmp1(1:ITMAXDcex,1:IDMAXDcex)=
     .   DCEX(ic,1:ITMAXDcex,1:IDMAXDcex)
        do cT=0,rtnt
          logT=log10(rtt(cT))
          do cD=0,rtnn
            logD=log10(rtn(cD)/1.0e6_R8)
            call int2d(DTMP1(1:itdimd,1:iddimd),itdimd,iddimd,
     1       DTEVDcex,ITMAXDcex,DDENSDcex,IDMAXDcex,
     2       logT,logD,ans)
            rtcx(cT,cD,ic)=10.0_R8**(ans-6.0_R8)
          enddo
        enddo
      enddo
c     filling rtzz(effective charge)
      write(*,*) 'rtzz'
      do ic=0,ncnct-1
        dtmp1(1:ITMAXDzav,1:IDMAXDzav)=
     .   DZAV(ic+1,1:ITMAXDzav,1:IDMAXDzav)
        do cT=0,rtnt
          logT=log10(rtt(cT))
          do cD=0,rtnn
            logD=log10(rtn(cD)/1.0e6_R8)
            call int2d(DTMP1(1:itdimd,1:iddimd),itdimd,iddimd,
     1       DTEVDzav,ITMAXDzav,DDENSDzav,IDMAXDzav,
     2       logT,logD,ans)
            rtzz(cT,cD,ic)=ans
          enddo
        enddo
      enddo
c     filling rtz2(effective square charge)
      write(*,*) 'rtz2'
      do ic=0,ncnct-1
        dtmp1(1:ITMAXDzav,1:IDMAXDzav)=
     .   DZA2(ic+1,1:ITMAXDzav,1:IDMAXDzav)
        do cT=0,rtnt
          logT=log10(rtt(cT))
          do cD=0,rtnn
            logD=log10(rtn(cD)/1.0e6_R8)
            call int2d(DTMP1(1:itdimd,1:iddimd),itdimd,iddimd,
     1       DTEVDza2,ITMAXDza2,DDENSDza2,IDMAXDza2,
     2       logT,logD,ans)
            rtz2(cT,cD,ic)=ans
          enddo
        enddo
      enddo
c     filling rtpt(cumulative ionisation potential), and rtpi(effective ionisation potential)
      write(*,*) 'rtpt,rtpi'
      do ic=0,ncnct-1
        do cT=0,rtnt
          dtmp1(1:ITMAXDpot,1:IDMAXDpot)=
     .     DPOT(ic+1,1:ITMAXDpot,1:IDMAXDpot)
          logT=log10(rtt(cT))
          do cD=0,rtnn
            logD=log10(rtn(cD)/1.0e6_R8)
            call int2d(DTMP1(1:itdimd,1:iddimd),itdimd,iddimd,
     1       DTEVDpot,ITMAXDpot,DDENSDpot,IDMAXDpot,
     2       logT,logD,ans)
            rtpi(cT,cD,ic)=ans
          enddo
        enddo
        if (ic.gt.0) then
          rtpt(0:rtnt,0:rtnn,ic)=
     1     rtpt(0:rtnt,0:rtnn,ic-1)+rtpi(0:rtnt,0:rtnn,ic)
        else
          rtpt(0:rtnt,0:rtnn,ic)=rtpi(0:rtnt,0:rtnn,ic)
        endif
      enddo
      
      call subend ()
      return

  990 continue
      call xerrab ('Check specifications in ratadas.filelist!')

      end subroutine ratadas

#ifdef ITM
      subroutine ratamns(nput,nc,rtt,rtnt,rtn,rtnn,
     &                   rtsa,rtra,rtqa,rtcx,rtrd,rtbr,
     &                   rtzz,rtz2,rtpt,rtpi,
     &                   iyear,fix_recomb,ncnct,izpmin,izpmax)
#ifdef GOT_AMNSPROTO
      use amns_types       ! IGNORE
!!      use amns_constants   ! IGNORE
      use amns_module      ! IGNORE
      use euitm_routines   ! IGNORE
#else
      use b2mod_types
#endif
      use b2mod_constants
     &     , only : qe
      use b2mod_solpstop

      implicit none

c     INPUT
      integer nput
      integer nc
      integer rtnt
      integer rtnn
      integer fix_recomb
      real (kind=R8) ::
     *   rtt(0:rtnt)
      real (kind=R8) ::
     *   rtn(0:rtnn)
c     OUTPUT:
      real (kind=R8) ::
     *   rtsa(0:rtnt,0:rtnn,0:nc)
      real (kind=R8) ::
     *   rtra(0:rtnt,0:rtnn,0:nc)
      real (kind=R8) ::
     *   rtqa(0:rtnt,0:rtnn,0:nc)
      real (kind=R8) ::
     *   rtcx(0:rtnt,0:rtnn,0:nc)
      real (kind=R8) ::
     *   rtrd(0:rtnt,0:rtnn,0:nc)
      real (kind=R8) ::
     *   rtbr(0:rtnt,0:rtnn,0:nc)
      real (kind=R8) ::
     *   rtzz(0:rtnt,0:rtnn,0:nc)
      real (kind=R8) ::
     *   rtz2(0:rtnt,0:rtnn,0:nc)
      real (kind=R8) ::
     *   rtpt(0:rtnt,0:rtnn,0:nc)
      real (kind=R8) ::
     *   rtpi(0:rtnt,0:rtnn,0:nc)
      integer iyear, ncnct
      integer izpmin(0:nc), izpmax(0:nc)
c     LOCAL VARIABLES:
      real (kind=R8) ::
     *   loc_te(0:rtnt,0:rtnn), loc_ne(0:rtnt,0:rtnn)
c     read in data
      integer ic, itnt, itnn


#ifdef GOT_AMNSPROTO
      type (amns_handle_type), save :: amns
      type (amns_handle_rx_type) ::
     &    amns_rc(0:nc), amns_ei(0:nc), amns_cx(0:nc), amns_br(0:nc),
     &    amns_lr(0:nc), amns_ze(0:nc), amns_ze2(0:nc), amns_eip(0:nc)
      type (amns_version_type) :: amns_database
      type (amns_reaction_type) :: rc_rx, ei_rx, cx_rx, br_rx,
     &     lr_rx, ze_rx, ze2_rx, eip_rx
      type (amns_reactants_type) :: species(0:nc)
      type (amns_query_type) :: query
      type (amns_answer_type) :: answer
      type (amns_set_type) :: set
      character (len=80) :: format
#endif

      logical silent_open
      external silent_open

*   ..subprogram start-up calls
      call subini ('ratamns')

#ifdef GOT_AMNSPROTO
      write(*,*) 'ITM AMNSPROTO data used (via UAL)'
      call ITM_AMNS_SETUP(amns)
      query%string='version'
      call ITM_AMNS_QUERY(amns,query,answer)
      write(*,*) 'AMNS data base version = ',trim(answer%string)
      rc_rx%string='RC'
      ei_rx%string='EI'
      cx_rx%string='CX'
      br_rx%string='BR'
      lr_rx%string='LR'
      ze_rx%string='ZE'
      ze2_rx%string='ZE2'
      eip_rx%string='EIP'
      format=
     & '(''ZN = '',f5.2,'', IS = '',i3,'', RX = '',a,'', SRC = '',a)'
      query%string='source'
      do ic=0, nc
         allocate(species(ic)%components(1))
         species(ic)%components(1)%ZN=nc
         species(ic)%components(1)%ZA=ic
         species(ic)%components(1)%MI=-1
      enddo
      do ic=0, nc-1
         call ITM_AMNS_SETUP_TABLE(amns, ei_rx,
     &        species(ic), amns_ei(ic))
         call ITM_AMNS_QUERY_TABLE(amns_ei(ic),query,answer)
         write(*,format) species(ic)%components(1)%ZN,
     &        ic, trim(ei_rx%string), trim(answer%string)
         call ITM_AMNS_SETUP_TABLE(amns, lr_rx,
     &        species(ic), amns_lr(ic))
         call ITM_AMNS_QUERY_TABLE(amns_lr(ic),query,answer)
         write(*,format) species(ic)%components(1)%ZN,
     &        ic, trim(lr_rx%string), trim(answer%string)
      enddo
      do ic=1, nc
         call ITM_AMNS_SETUP_TABLE(amns, rc_rx,
     &        species(ic), amns_rc(ic))
         call ITM_AMNS_QUERY_TABLE(amns_rc(ic),query,answer)
         write(*,format) species(ic)%components(1)%ZN,
     &        ic, trim(rc_rx%string), trim(answer%string)
         call ITM_AMNS_SETUP_TABLE(amns, br_rx,
     &        species(ic), amns_br(ic))
         call ITM_AMNS_QUERY_TABLE(amns_br(ic),query,answer)
         write(*,format) species(ic)%components(1)%ZN,
     &        ic, trim(br_rx%string), trim(answer%string)
         call ITM_AMNS_SETUP_TABLE(amns, cx_rx,
     &        species(ic), amns_cx(ic))
         call ITM_AMNS_QUERY_TABLE(amns_cx(ic),query,answer)
         write(*,format) species(ic)%components(1)%ZN,
     &        ic, trim(cx_rx%string), trim(answer%string)
      enddo
      do ic=0, nc
         call ITM_AMNS_SETUP_TABLE(amns, ze_rx,
     &        species(ic), amns_ze(ic))
         call ITM_AMNS_QUERY_TABLE(amns_ze(ic),query,answer)
         write(*,format) species(ic)%components(1)%ZN,
     &        ic, trim(ze_rx%string), trim(answer%string)
         call ITM_AMNS_SETUP_TABLE(amns, ze2_rx,
     &        species(ic), amns_ze2(ic))
         call ITM_AMNS_QUERY_TABLE(amns_ze2(ic),query,answer)
         write(*,format) species(ic)%components(1)%ZN,
     &        ic, trim(ze2_rx%string), trim(answer%string)
         call ITM_AMNS_SETUP_TABLE(amns, eip_rx,
     &        species(ic), amns_eip(ic))
         call ITM_AMNS_QUERY_TABLE(amns_eip(ic),query,answer)
         write(*,format) species(ic)%components(1)%ZN,
     &        ic, trim(eip_rx%string), trim(answer%string)
      enddo
!!!    stop
#else
      write(*,*) 'No ITM-AMNS data available'
#endif
      do itnn=0, rtnn
         do itnt=0, rtnt
            loc_te(itnt, itnn) = rtt(itnt)
            loc_ne(itnt, itnn) = rtn(itnn)
         enddo
      enddo
      do ic = 0, nc
         izpmin(ic) = ic
         izpmax(ic) = ic
      enddo
      ncnct=nc+1

      rtsa=0.0_R8
      rtra=0.0_R8
      rtqa=0.0_R8
      rtcx=0.0_R8
      rtrd=0.0_R8
      rtbr=0.0_R8
      rtzz=-1.0_R8
      rtz2=-1.0_R8
      rtpi=-1.0_R8
      if(fix_recomb.ne.0) then
        write(*,*) 'RATAMNS: FIX_RECOMB <> 0 : ',fix_recomb
      endif

#ifdef GOT_AMNSPROTO
c     filling rtsa(ionisation)
      write(*,*) 'rtsa'
      do ic=0,nc-1
         call ITM_AMNS_RX(amns_ei(ic),rtsa(:,:,ic),loc_te,loc_ne)
         call ITM_AMNS_FINISH_TABLE(amns_ei(ic))
      enddo
c     filling rtra(recombination)
      write(*,*) 'rtra'
      do ic=1,nc
         call ITM_AMNS_RX(amns_rc(ic),rtra(:,:,ic),loc_te,loc_ne)
         call ITM_AMNS_FINISH_TABLE(amns_rc(ic))
      enddo
c     filling rtcx(charge exchange)
      write(*,*) 'rtcx'
      do ic=1,nc
         call ITM_AMNS_RX(amns_cx(ic),rtcx(:,:,ic),loc_te,loc_ne)
         call ITM_AMNS_FINISH_TABLE(amns_cx(ic))
      enddo
c     filling rtzz(effective charge)
      write(*,*) 'rtzz'
      do ic=0,nc
         call ITM_AMNS_RX(amns_ze(ic),rtzz(:,:,ic),loc_te,loc_ne)
         call ITM_AMNS_FINISH_TABLE(amns_ze(ic))
      enddo
c     filling rtz2(effective square charge)
      write(*,*) 'rtz2'
      do ic=0,nc
         call ITM_AMNS_RX(amns_ze2(ic),rtz2(:,:,ic),loc_te,loc_ne)
         call ITM_AMNS_FINISH_TABLE(amns_ze2(ic))
      enddo
c     filling rtpt(cumulative ionisation potential), and rtpi(effective ionisation potential)
      write(*,*) 'rtpt,rtpi'
      do ic=0,nc
         call ITM_AMNS_RX(amns_eip(ic),rtpi(:,:,ic),loc_te,loc_ne)
         call ITM_AMNS_FINISH_TABLE(amns_eip(ic))
         if (ic.gt.0) then
            rtpt(0:rtnt,0:rtnn,ic)=
     1           rtpt(0:rtnt,0:rtnn,ic-1)+rtpi(0:rtnt,0:rtnn,ic)
         else
            rtpt(0:rtnt,0:rtnn,ic)=rtpi(0:rtnt,0:rtnn,ic)
         endif
      enddo
c     filling rtqa(electron heat loss), rtrd(line radiation) and rtbr(bremsstrahlung)
      write(*,*) 'rtqa,rtrd,rtbr'
      do ic=0,nc
         do itnn=0, rtnn
            do itnt=0, rtnt
               rtbr(itnt,itnn,ic) =
     &              1.69e-38_R8/qe*ic**2*sqrt(rtt(itnt))
            enddo
         enddo
      enddo
      do ic=0,nc-1
         call ITM_AMNS_RX(amns_lr(ic),rtrd(:,:,ic),loc_te,loc_ne)
         call ITM_AMNS_FINISH_TABLE(amns_lr(ic))
      enddo
      rtrd=rtrd/qe
      if(fix_recomb.eq.0) then
         rtqa(:,:,:)=rtbr(:,:,:)
      else
         call xerrab ('RATAMNS: fix_recomb <> 0 : Not yet coded!')
      endif
      rtqa(:,:,:)=rtqa(:,:,:)+rtrd(:,:,:)
      do ic=0,nc-1
         rtqa(:,:,ic)=rtqa(:,:,ic)+rtsa(:,:,ic)*rtpi(:,:,ic+1)
      enddo
      call ITM_AMNS_FINISH(amns)
#endif

      call subend ()
      return
      end subroutine ratamns
#endif
