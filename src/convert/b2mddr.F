      subroutine b2mddr (ninp, nout, ns, mpg, geo, state, state_ext)
      use b2mod_types
      use b2mod_version
      use b2mod_math
      use b2mod_ppout
      use b2mod_mdsplus
      use b2mod_b2cmfs
      use b2mod_b2cmrc
      use b2mod_b2cmgs
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_b2cmpt
      use b2mod_b2cmwg
      use b2mod_subsys
      use b2mod_ipmain
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2us_io
      use b2mod_switches
      implicit none
*   ..input arguments (unchanged on exit)
      integer ninp(0:6), nout(0:2), ns
*   ..main data structures
      type(mapping)    :: mpg
      type(geometry)   :: geo
      type(B2State)    :: state
      type(B2StateExt) :: state_ext
      type(switches)   :: switch
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
#ifdef MDSPLUS
#     include "mdslib.inc"
      character :: tcl*80, shottmp*12, caseline*32
      integer :: status, shotnumber, status2, pos_index,
     1 in_jxi, in_jxa, in_jsep, retlen, svn_revision
      character :: mds_server*(128)
      character*256 filename
      logical :: file_exists
#endif
*   ..local variables
      integer k, is, idum(0:9), ismain
      integer external_species
      character lblgm*120, lblcp*120, lblmn*120, lblrc*120
      character cnamip*80, cvalip*80, line*256, ligne*256
*   ..initialisation
      save ismain, external_species
      data ismain/1/, external_species/0/
*   ..procedures
      logical is_comment
      external is_comment
      external xertst, xerrab, cfruch, cfruin
      external b2rups, b2rugm, b2rucp, b2mds, ipgeti
      external b2rurc, b2spel, b2sqel, b2sqcx, b2spcx, b2xvps
      external b2xxid, b2rflb, b2rfcp
*   ..namelist
      character :: exp*32, comment*128
      integer :: shot, overwrite_shotnumber
      real (kind=R8) :: time
      logical :: timedep, snapshot, tallies, movies
      namelist /b2md_namelist/ exp, shot, time, comment,
     1 timedep, snapshot, tallies, movies, overwrite_shotnumber

*   ..subprogram start-up calls
      call subini ('b2mddr')
*   ..test ninp, nout
      call xertst (1.le.ninp(0).and.1.le.ninp(1).and.1.le.ninp(2).and.
     &             1.le.ninp(3).and.1.le.ninp(5).and.1.le.ninp(6).and.
     &             1.le.nout(0).and.1.le.nout(1).and.1.le.nout(2),
     &  'faulty argument ninp, nout')
*   ..test dimensions
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (ns.le.nsdecl, 'faulty parameter nsdecl')

* ..computation
*   ..obtain lblmn
      call b2xxid ('b2mn', lblmn(1:60))
      call b2rflb (ninp(5), lblmn(61:120))
*   ..read revised parameters in physics common
*     (ninp(0) need supply only parameters that are to be changed from
*     the values read in as default physics common.)
      call b2rfcp (ninp(5), ns)
*   ..read and echo code internal parameters
      write (nout(0),'(/2x,a)') 'non-default internal parameters:'
    1 continue
       read (ninp(5),'(a80)',end=2,err=91) line
       if (.not.is_comment(line)) then
         ligne = line
         call strip_spaces(ligne)
         if (ligne(1:1).eq.'''') then
           read (line,*) cnamip, cvalip
           write (nout(0),'(4x,a,2x,a)') cnamip, trim(cvalip)
           call ipsetc (cnamip, cvalip)
         endif
       endif
       goto 1
    2 continue
      write (nout(0),'(2x,a)') '(end of list of internal parameters)'
      call read_first_switches(switch)
      call set_defaults_switches(switch,ns)
      call read_switches(switch)
      call check_switches(switch,ns)
*   ..read geometry
      call cfruch (ninp(1), 120, lblgm, 'label')
      call cfruin (ninp(1), 1, idum, 'isymm')
      isymm = idum(0)
      call read_b2fgmtry (ninp(1), mpg, geo)
      call init_geometry (mpg, geo, switch)
*   ..read physics common
      call cfruch (ninp(2), 120, lblcp, 'label')
      call b2rucp (ninp(2), b2fparam_version, ns)
*   ..read plasma state file label
      call cfruch (ninp(3), 120, lblmn, 'label')
*   ..read species data
      call b2ruzd (ninp(3), b2fstate_version,
     &  ns, zamin, zamax, zn, am, .true.)
*   ..read plasma state
      call read_b2fstate (ninp(3), mpg%nCv, mpg%nFc, ns, state)
*   ..test plasma state and fluxes
      call b2xvps (mpg%nCv, mpg%nFc, ns, state%pl, state%dv)
*   ..read atomic rate data
      call cfruch (ninp(6), 120, lblrc, 'label')
      adpak_used = index(lblrc,'ADPAK').gt.0
      call read_b2mod_ppout('rates', ninp(6), lblrc)
      call b2rurc (ninp(6), b2frates_version)
*   ..test rtns, rtzmin, rtzmax, rtzn
      call xertst (ns.eq.rtns, 'faulty input rtns--wrong rate table?')
      do is = 0, ns-1
       call xertst (ltst(zamin(is),rtzmin(is)),
     &   'faulty input rtzmin--wrong rate table?')
       call xertst (ltst(zamax(is),rtzmax(is)),
     &   'faulty input rtzmax--wrong rate table?')
       call xertst (ltst(zn(is),rtzn(is)),
     &   'faulty input rtzn--wrong rate table?')
      enddo
*   ..compute log-log electron rate coefficients
      call b2spel (mpg%nCv, ns, ev, state%pl%te, state%dv%ne, state%rt)
*   ..compute electron rate coefficients
      call ipgeti ('b2mndr_ismain', ismain)
      call xertst (0.le.ismain.and.ismain.lt.ns,
     . 'invalid main plasma species index ismain')
      call xertst (.not.is_neutral(ismain),
     . 'invalid main plasma species ismain; must not be neutral')
      call b2sqel (mpg%nCv, ns, ismain, switch, ev,
     .  state%pl%te, state%rt, state%rtw)
*   ..find nscx, iscx
*     (indices for neutral hydrogen species)
      state%rt%nscx = 0
      do is = 0, ns-1
        if (is_neutral(is).and.nint(zn(is)).eq.1) then
          call xertst (state%rt%nscx.lt.state%rt%nscxmax,
     1     'too many neutral hydrogenic species')
          state%rt%iscx(state%rt%nscx) = is
          state%rt%nscx = state%rt%nscx+1
         endif
       enddo
       do k = state%rt%nscx, state%rt%nscxmax-1
         state%rt%iscx(k) = -1
       enddo
*    ..compute log-log charge exchange rate coefficients
      do k = 0, state%rt%nscx-1
        call b2spcx (mpg%nCv, ns, ev, am(state%rt%iscx(k)),
     &   state%pl%ti, state%dv%ne,
     &   state%rt%rlcx(1:mpg%nCv,0:1,0:ns-1,k))
      enddo
*   ..compute charge exchange rate coefficients
      do k = 0, state%rt%nscx-1
        call b2sqcx (mpg%nCv, ns, ev, switch, am(state%rt%iscx(k)),
     &   state%pl%ti, state%pl%tn,
     &   state%rt%rlcx(1:mpg%nCv,0:1,0:ns-1,k), state%rtw%rcx)
      enddo
*   ..obtain external species data if any
      call ipgeti ('b2mndr_external_species', external_species)
#ifdef WG_TODO
      if (external_species.ne.0) then
        call read_b2mod_external(nx,ny,ns,ismain)
      endif
#endif

*   ..process
      exp='NOT_SET'
      comment='NOT_SET'
      shot=0
      time=0.0_R8
      timedep=.true.
      snapshot=.true.
      tallies=.true.
      movies=.true.
      overwrite_shotnumber=0
      read(ninp(0),b2md_namelist,err=99)
      write(nout(0),b2md_namelist)

#ifdef MDSPLUS
      read(*,'(a120)') mds_server
      if(mds_server(1:10).ne.'MDSSERVER:') then
        write(*,*) 'Invalid server label from "mds_id"'
        call xerrab ('Invalid server label')
      else
        mds_server(1:128)=mds_server(11:128)
      endif

      write(nout(0),'(a,a)')
     1 'Test implementation of writing the SOLPS output ',
     2 'to a MDSPLUS server'
      write(nout(0),'(a,a)') 'The server used was ',
     1 trim(mds_server)

      status=mdsconnect(trim(mds_server)//char(0))
      write(*,*) 'mdsconnect'
      if(status.lt.0) then
        write(*,*) ' return status: ',status
        call xertst (status.ge.0,
     .   'Error connecting to the MDSplus server ' // trim(mds_server))
      endif

      if(overwrite_shotnumber.le.0) then

        status=mdsvalue2(
     *   'CURRENT_SHOT("SOLPS")'//char(0),
     *   descr2(IDTYPE_LONG,0),shotnumber,
     *   0,retlen)
        if(mod(status,2).eq.1) then
          write(*,*) 'old shotnumber: ',shotnumber
          write(nout(0),'(a,i8)') 'The previous shot was # ',shotnumber
          call xertst (shotnumber.ne.0,
     *     'Stopping because previous shotnumber = 0')
        else
          call mds_err_msg(status)
          call xerrab ('Invalid status value')
        endif

        status=mdsvalue2(
     *   'TCL("set current solps/increment"),CURRENT_SHOT("SOLPS")'
     *   //char(0),
     *   descr2(IDTYPE_LONG,0),shotnumber,
     *   0,retlen)
        if(mod(status,2).eq.1) then
          call mds_err_msg(status)
          write(*,*) 'new shotnumber: ',shotnumber
          write(nout(0),'(a,i8)') 'The new shot was # ',shotnumber
          call xertst (shotnumber.ne.0,
     *     'Stopping because shotnumber = 0')
          dbsolpsindex=shotnumber
        else
          call mds_err_msg(status)
          call xerrab ('Invalid status value')
        endif

        tcl='TCL("set tree solps")'
        status=mdsvalue2(
     *   trim(tcl)//char(0),
     *   descr2(IDTYPE_LONG,0),status2,
     *   0,retlen)
        call mds_err_msg(status)
        call mds_err_msg(status2)

        write(shottmp,'(i12)') shotnumber
        tcl='TCL("create pulse ' // shottmp // '")'
        status=mdsvalue2(
     *   trim(tcl)//char(0),
     *   descr2(IDTYPE_LONG,0),status2,
     *   0,retlen)
        call mds_err_msg(status)
        call mds_err_msg(status2)

cdpc        status=mdsvalue2(
cdpc     *   'CURRENT_SHOT("SOLPS")'//char(0),
cdpc     *   descr2(IDTYPE_LONG,0),shotnumber,
cdpc     *   0,retlen)
cdpc        if(mod(status,2).eq.1) then
cdpc          write(*,*) 'current shotnumber: ',shotnumber
cdpc        else
cdpc          call mds_err_msg(status)
cdpc          call xerrab ('Invalid status value')
cdpc        endif

      else

        shotnumber=overwrite_shotnumber
        write(*,*) 'Overwriting: ',shotnumber
        write(nout(0),'(a,i8)') 'Overwriting shot # ',shotnumber
        call xertst (shotnumber.ne.0,
     *   'Stopping because overwriting shotnumber = 0')
        dbsolpsindex=shotnumber

        tcl='TCL("set tree solps")'
        status=mdsvalue2(
     *   trim(tcl)//char(0),
     *   descr2(IDTYPE_LONG,0),status2,
     *   0,retlen)
        call mds_err_msg(status)
        call mds_err_msg(status2)

        write(shottmp,'(i12)') shotnumber
        tcl='TCL("create pulse ' // shottmp // '")'
        status=mdsvalue2(
     *   trim(tcl)//char(0),
     *   descr2(IDTYPE_LONG,0),status2,
     *   0,retlen)
        call mds_err_msg(status)
        call mds_err_msg(status2)

      endif

      filename='shotnumber.history'
      call find_file(filename,file_exists)
      open(99,file=filename,position='append')
      write(99,'(i8)') shotnumber
      close(99)

      status=mdsopen('solps'//char(0),shotnumber)
      write(*,*) 'mdsopen'
      call mds_err_msg(status)

   10 continue
      read(*,'(a256)',end=20) line
      pos_index=index(line,':')
      if(pos_index.eq.0) then
        write(*,*) 'Line "',trim(line),'" not recognised'
      else
        caseline=line(1:pos_index-1)
        line(1:256)=line(pos_index+1:256)
        select case (caseline)
          case('date')
            if(len_trim(line).eq.0) line='UNKNOWN DATE' !
            call mds_put_string('\TOP.IDENT:DATE',line)
            write(nout(0),'(a,a)') '         DATE : ',trim(line)
            dbdate=trim(line)
          case('dir')
            if(len_trim(line).eq.0) line='UNKNOWN DIRECTORY'
            call mds_put_string('\TOP.IDENT:DIRECTORY',line)
            write(nout(0),'(a,a)') '    DIRECTORY : ',trim(line)
            dbdirectory=trim(line)
          case('user')
            if(len_trim(line).eq.0) line='UNKNOWN USER'
            call mds_put_string('\TOP.IDENT:USER',line)
            write(nout(0),'(a,a)') '         USER : ',trim(line)
            dbuser=trim(line)
          case('hostname')
            if(len_trim(line).eq.0) line='UNKNOWN HOST'
            call mds_put_string('\TOP.IDENT:HOST',line)
            write(nout(0),'(a,a)') '         HOST : ',trim(line)
            dbhost=trim(line)
          case('solpstop')
            if(len_trim(line).eq.0) line='UNKNOWN SOLPSTOP'
            call mds_put_string('\TOP.IDENT:SOLPSTOP',line)
            write(nout(0),'(a,a)') '     SOLPSTOP : ',trim(line)
            dbsolpstop=trim(line)
          case('objectcode')
            if(len_trim(line).eq.0) line='UNKNOWN OBJECTCODE'
            call mds_put_string('\TOP.IDENT:OBJECTCODE',line)
            write(nout(0),'(a,a)') '   OBJECTCODE : ',trim(line)
            dbobjectcode=trim(line)
          case('b2fstatiid')
            if(len_trim(line).eq.0) then
              line='UNKNOWN MD5SUM'
            elseif(line(34:41).ne.'b2fstati') then
              write(*,'(3a)') '"',line(34:41),'"'
              line='MISMATCHED MD5SUM'
            else
              line=line(1:32)
            endif
            call mds_put_string('\TOP.IDENT:B2FSTATIID',line)
            write(nout(0),'(a,a)') '   B2FSTATIID : ',trim(line)
            dbb2fstatiid=trim(line)
          case('b2fstateid')
            if(len_trim(line).eq.0) then
              line='UNKNOWN MD5SUM'
            elseif(line(34:41).ne.'b2fstate') then
              write(*,'(3a)') '"',line(34:41),'"'
              line='MISMATCHED MD5SUM'
            else
              line=line(1:32)
            endif
            call mds_put_string('\TOP.IDENT:B2FSTATEID',line)
            write(nout(0),'(a,a)') '   B2FSTATEID : ',trim(line)
            dbb2fstateid=trim(line)
          case('b2fplasmaid')
            if(len_trim(line).eq.0) then
              line='UNKNOWN MD5SUM'
            elseif(line(34:45).ne.'../b2fplasma') then
              write(*,'(3a)') '"',line(34:45),'"'
              line='MISMATCHED MD5SUM'
            else
              line=line(1:32)
            endif
            call mds_put_string('\TOP.IDENT:B2FPLASMAID',line)
            write(nout(0),'(a,a)') '  B2FPLASMAID : ',trim(line)
            dbb2fplasmaid=trim(line)
          case('b2fgmtryid')
            if(len_trim(line).eq.0) then
              line='UNKNOWN MD5SUM'
            elseif(line(34:41).ne.'b2fgmtry') then
              write(*,'(3a)') '"',line(34:41),'"'
              line='MISMATCHED MD5SUM'
            else
              line=line(1:32)
            endif
            call mds_put_string('\TOP.IDENT:B2FGMTRYID',line)
            write(nout(0),'(a,a)') '   B2FGMTRYID : ',trim(line)
            dbb2fgmtryid=trim(line)
          case('b2fratesid')
            if(len_trim(line).eq.0) then
              line='UNKNOWN MD5SUM'
            elseif(line(34:41).ne.'b2frates') then
              write(*,'(3a)') '"',line(34:41),'"'
              line='MISMATCHED MD5SUM'
            else
              line=line(1:32)
            endif
            call mds_put_string('\TOP.IDENT:B2FRATESID',line)
            write(nout(0),'(a,a)') '   B2FRATESID : ',trim(line)
            dbb2fratesid=trim(line)
          case('jxijxajsep')
            in_jxi=-2
            in_jxa=-2
            in_jsep=-2
            read(line,*,err=10,end=20) in_jxi,in_jxa,in_jsep
          case('Revision')
            read(line,*,err=10,end=20) svn_revision
            call mds_put_int('\TOP.IDENT.VERSION:SVN',svn_revision)
            write(nout(0),'(a,i8)') '          SVN : ',svn_revision
            dbsvn=svn_revision
          case('git')
            call mds_put_string('\TOP.IDENT.VERSION:GIT',trim(line))
            write(nout(0),'(a,a)') '          GIT : ',trim(line)
            dbgit=trim(line)
          case default
            write(*,*) 'Line option "',trim(caseline),'" not recognised'
        end select
      endif
      goto 10

   20 continue

      call mds_put_string('\TOP.IDENT:GEOMETRY',lblgm)
      write(nout(0),'(a,a)') '     GEOMETRY : ',trim(lblgm)
      dbgeometry=trim(lblgm)
      call mds_put_string('\TOP.IDENT:PHYSICS',lblcp)
      write(nout(0),'(a,a)') '      PHYSICS : ',trim(lblcp)
      dbphysics=trim(lblcp)
      call mds_put_string('\TOP.IDENT:MAIN',lblmn)
      write(nout(0),'(a,a)') '         MAIN : ',trim(lblmn)
      dbmain=trim(lblmn)

      call mds_put_string('\TOP.IDENT:COMMENTS',comment)
      dbcomments=comment
      call mds_put_string('\TOP.IDENT:EXP',exp)
      dbexp=exp
      call mds_put_int('\TOP.IDENT:SHOT',shot)
      dbshot=shot
      call mds_put_double('\TOP.IDENT:TIME','s',time)
      dbtime=time

      call mds_put_int('\TOP.IDENT:NCV',mpg%nCv)
      dbnCv=mpg%nCv
      call mds_put_int('\TOP.IDENT:NFC',mpg%nFc)
      dbnFc=mpg%nFc
      call mds_put_int('\TOP.IDENT:NS',ns)
      dbns=ns

      call mds_put_string('\TOP.IDENT.VERSION:b2fstati',
     *     b2fstati_version)
      write(nout(0),'(a,a)') '   V_B2FSTATI : ',b2fstati_version
      db_version_b2fstati=b2fstati_version
      call mds_put_string('\TOP.IDENT.VERSION:b2fstate',
     *     b2fstate_version)
      write(nout(0),'(a,a)') '   V_B2FSTATE : ',b2fstate_version
      db_version_b2fstate=b2fstate_version
      call mds_put_string('\TOP.IDENT.VERSION:b2frates',
     *     b2frates_version)
      write(nout(0),'(a,a)') '   V_B2FRATES : ',b2frates_version
      db_version_b2frates=b2frates_version
      call mds_put_string('\TOP.IDENT.VERSION:b2fmovie',
     *     b2fmovie_version)
      write(nout(0),'(a,a)') '   V_B2FMOVIE : ',b2fmovie_version
      db_version_b2fmovie=b2fmovie_version
      call mds_put_string('\TOP.IDENT.VERSION:b2ftrack',
     *     b2ftrack_version)
      write(nout(0),'(a,a)') '   V_B2FTRACK : ',b2ftrack_version
      db_version_b2ftrack=b2ftrack_version
      call mds_put_string('\TOP.IDENT.VERSION:plasmastate',
     *     plasmastate_version)
      write(nout(0),'(a,a)') 'V_PLASMASTATE : ',plasmastate_version
      db_version_plasmastate=plasmastate_version
      call mds_put_string('\TOP.IDENT.VERSION:b2fplasma',
     *     b2fplasma_version)
      write(nout(0),'(a,a)') '  V_B2FPLASMA : ',b2fplasma_version
      db_version_b2fplasma=b2fplasma_version

      if(timedep) then
        filename='b2time.nc'
        call find_file(filename,file_exists)
        if(file_exists) then
          write(*,*) 'Calling b2cdf2mds'
          call b2cdf2mds
          write(*,*) 'Returned from b2cdf2mds'
        else
          write(*,*) 'Skipped calling b2cdf2mds because '//
     .               'no b2time.nc file found'
        endif
      else
        write(*,*) 'Skipped calling b2cdf2mds because timedep=F'
      endif
#ifdef WG_TODO
      if(snapshot) then
        write(*,*) 'Calling b2mds'
        call b2mds(nx,ny,ns,in_jxi,in_jxa,in_jsep,switch%BoRiS)
        write(*,*) 'Returned from b2mds'
      else
        write(*,*) 'Skipped calling b2mds because snapshot=F'
      endif
#endif
      if(tallies) then
        filename='b2tallies.nc'
        call find_file(filename,file_exists)
        if(file_exists) then
          write(*,*) 'Calling b2cdftallies2mds'
          call b2cdftallies2mds
          write(*,*) 'Returned from b2cdftallies2mds'
        else
          write(*,*)
     1     'Skipped calling b2cdftallies2mds because '//
     2     'no b2tallies.nc file found'
        endif
      else
        write(*,*) 'Skipped calling b2cdftallies2mds because tallies=F'
      endif
      if(movies) then
        filename='b2movies.nc'
        call find_file(filename,file_exists)
        if(file_exists) then
          write(*,*) 'Calling b2cdfmovies2mds'
          call b2cdfmovies2mds
          write(*,*) 'Returned from b2cdfmovies2mds'
        else
          write(*,*)
     1     'Skipped calling b2cdfmovies2mds because '//
     2     'no b2movies.nc file found'
        endif
      else
        write(*,*) 'Skipped calling b2cdfmovies2mds because movies=F'
      endif

      status=mdsclose('solps'//char(0),shotnumber)
      write(*,*) 'mdsclose'
      call mds_err_msg(status)
      call b2mod_mdsplus_write_to_db(nout(0))
      status=mdsdisconnect()
      write(*,*) 'mdsdisconnect'

* ..return
#else
      call xerrab('b2mddr: ' //
     1 'MDSPLUS routine called when not compiled with MDSPLUS')
#endif
      call dealloc_b2mod_b2cmrc
      call subend ()
      return

*-----------------------------------------------------------------------
*.error conditions

   91 call xerrab ('error trying to read internal parameters')
   99 call xerrab ('error trying to read b2md_namelist')

*-----------------------------------------------------------------------
*.end b2mddr

      end subroutine b2mddr

