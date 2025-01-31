      subroutine b2paspwr(ic)
c
c  version : 13.08.98 13:36
c
c======================================================================
c***  Adds a line to the file ../sum_tmp/spec_pwr.dat with the power
c***  entering both divertors, power radiated along the divertor legs,
c***  power radiated near the targets, areas of the separatrix branches
c***  and mean specific power radiated along the divertor legs per
c***  square meter. ig-1000 cells are treated as "close to the targets"
c======================================================================
      use helper
      use b2mod_types
      use b2mod_diag
      use b2mod_ipmain
      use b2mod_b2plot
      use b2mod_b2plot_wall_loading
      use b2us_data
      implicit none
c
      integer, intent(in) :: ic
      integer lun,m_hlp,m_hlpa
      parameter (lun=96, m_hlp=9, m_hlpa=1)
      integer, save :: shot = 0, run = 0
      character*12 case_number
      character*16 username
      character*24 database,version
      character*256 filename,path,imasdir,home_dir,zeile
      integer i,j,k,l,m,lv_hlp,n_hlp
      integer iCb,iCn,iCp,iCv,iCx,iFc,iFn
      real (kind=R8) ::
     ,  hlp_pwr_in(ntrgsx),hlp_rad_div(ntrgsx),
     ,  hlp_rad_trg(ntrgsx),hlp_area(ntrgsx),
     ,  hlp_spwr(ntrgsx),hlp_spwt(ntrgsx),
     ,  hlp_rdfr(ntrgsx),hlp_rad_tot(ntrgsx),
     ,  hlp_area_t(ntrgsx)
      real (kind=R8) ::
     ,  equ_hlp(ntrgsx,m_hlp),hlp_r,u,eqa_hlp(m_hlpa),hlp_tot
      equivalence
     ,  (hlp_tot,eqa_hlp(1))
      equivalence
     ,  (hlp_pwr_in,equ_hlp(1,1)),  (hlp_rad_div,equ_hlp(1,2)),
     ,  (hlp_rad_trg,equ_hlp(1,3)), (hlp_area,equ_hlp(1,4)),
     ,  (hlp_spwr,equ_hlp(1,5)),    (hlp_spwt,equ_hlp(1,6)),
     ,  (hlp_rdfr,equ_hlp(1,7)),    (hlp_rad_tot,equ_hlp(1,8)),
     ,  (hlp_area_t,equ_hlp(1,9))
      logical ex,not_default
      character hlp_fn*24,hlp_c(ntrgsx)*2,hlp_u,hlp_frm*256,hlp_frh*256,
     ,  hlp_hdr(m_hlp)*4,hlp_hda(m_hlpa)*4,hlp_inp*8,hlp_s(ntrgsx)
      data hlp_fn / 'spec_pwr.dat' /
      data hlp_c / 'i ', 'iu', 'ou', 'o ' /
      data hlp_hdr / 'Pin_','Pdv_','Ptg_','Sdv_','Psw_',
     ,               'Pst_','RdF_','Prd_','Sdt_' /
      data hlp_hda / 'Pin' /
      logical streql
      character*16 usrnam
      character*24 device_env
      external streql,usrnam
#ifndef NO_GETENV
      integer lenval,ierror
#ifndef USE_PXFGETENV
      intrinsic get_environment_variable
#endif
#endif
c======================================================================

      if (ntrgts.eq.0) ntrgts=maxval(mpg%strDiv)
      n_hlp=ntrgts*m_hlp+m_hlpa

c*** Calculate the input power for each divertor
c*** Calculate the areas and power radiated in the divertors
c*** Calculate the power radiated near the targets

      hlp_tot=0.0_R8
      do l=1,ntrgts
        m=0
        hlp_area(l)=0.0_R8
        hlp_area_t(l)=0.0_R8
        hlp_pwr_in(l)=0.0_R8
        hlp_rad_div(l)=0.0_R8
        hlp_rad_trg(l)=0.0_R8
        do j=mpg%ifdiv(l),mpg%divFcP(l,2)
          iFc = mpg%divFc(mpg%divFcP(l,1) + j - 1)
          iCv = min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
          iCn=iCv
          iCb=iCv
          iCx=0
          do while (iCx.eq.0)
            iCp = 0
            m = m+1
            do i = mpg%cvFcP(iCn,1),
     .             mpg%cvFcP(iCn,1) + mpg%cvFcP(iCn,2) - 1
              iFc = mpg%cvFc(i)
              if (j.eq.mpg%ifdiv(l)) then
                if (mpg%fcFs(iFc).eq.mpg%iFssep .or.
     .             (mpg%fcFs(iFc).eq.mpg%iFssep2.and.mpg%iFssep2.ne.0))
     >           then
                  if (m.ge.ic) then
                    hlp_area(l)=hlp_area(l)+geo%fcS(iFc)
                  else
                    hlp_area_t(l)=hlp_area_t(l)+geo%fcS(iFc)
                  end if
                end if
              end if
              if (iCp.ne.0) cycle
              if (mpg%fcFs(iFc).ne.0) cycle
              if (mpg%fcCv(iFc,1).gt.mpg%nCi.or.
     .            mpg%fcCv(iFc,2).gt.mpg%nCi) cycle
              if (mpg%fcCv(iFc,1).eq.iCn .and. mpg%fcCv(iFc,2).ne.iCb)
     >         then
                u = -1.0_R8
                iFn = iFc
                iCp = mpg%fcCv(iFn,2)
              end if
              if (mpg%fcCv(iFc,2).eq.iCn .and. mpg%fcCv(iFc,1).ne.iCb)
     >         then
                u = 1.0_R8
                iFn = iFc
                iCp = mpg%fcCv(iFn,1)
              end if
            end do
            iCb = iCn
            iCn = iCp
            do k = 0, state%ns-1
              if (m.ge.ic) then
                hlp_rad_div(l)=hlp_rad_div(l)+state%srw%rqbrm(iCb,k)
                hlp_rad_div(l)=hlp_rad_div(l)+state%srw%rqrad(iCb,k)
              else
                hlp_rad_trg(l)=hlp_rad_trg(l)+state%srw%rqbrm(iCb,k)
                hlp_rad_trg(l)=hlp_rad_trg(l)+state%srw%rqrad(iCb,k)
              end if
            end do
            if (mpg%cvReg(iCb).ne.mpg%cvReg(iCn)) then
              iCx = iCn
              hlp_pwr_in(l)=hlp_pwr_in(l)+
     .         u*(state%dv%fhe(iFn,0)+state%dv%fhi(iFn,0))
            end if
          end do
        end do
        hlp_tot=hlp_tot+hlp_pwr_in(l)
        hlp_pwr_in(l)=1.0e-6_R8*hlp_pwr_in(l)
        hlp_rad_div(l)=1.0e-6_R8*hlp_rad_div(l)
        hlp_rad_trg(l)=1.0e-6_R8*hlp_rad_trg(l)
      end do
      hlp_tot=1.e-6*hlp_tot
c      print *,'Areas: ',hlp_area,hlp_area_t
c      print *,'Input powers:',hlp_tot,hlp_pwr_in
c      print *,'Div. rad.: ',hlp_rad_div
c      print *,'Trg. rad.: ',hlp_rad_trg

c*** Total radiated power and its fraction

      do i=1,ntrgts
        hlp_rad_tot(i)=hlp_rad_div(i)+hlp_rad_trg(i)
        if(hlp_pwr_in(i).gt.0.0_R8) then
          hlp_rdfr(i)=hlp_rad_tot(i)/hlp_pwr_in(i)
        else
          hlp_rdfr(i)=0.0_R8
        end if
      end do

c*** Calculate the specific powers

      do l=1,ntrgts
        if(hlp_area(l).gt.0.0_R8) then
          hlp_spwr(l)=hlp_rad_div(l)/hlp_area(l)
          hlp_spwt(l)=hlp_rad_tot(l)/hlp_area_t(l)
        else
          hlp_spwr(l)=0.0_R8
          hlp_spwt(l)=0.0_R8
        end if
      end do
c      print *,'spwr : ',hlp_spwr

c*** Produce the output

      call find_file(hlp_fn,ex)
      open(lun,file=hlp_fn,err=900)
      write(hlp_frm,'(a14,i3,a9)') '(1x,a12,1x,1p,',n_hlp,'e12.4,i5)'
c      print *,'hlp_frm: ',trim(hlp_frm)
      if(ex) then

c*** read the file

        l=0
   10   read(lun,'(a1)',end=20) hlp_u
        if(hlp_u.eq.'#') then
          l=l+1
          go to 10
        end if
   20   if(l.eq.0) then
          write(0,*) 'b2paspwr: wrong format of the existing file ',
     ,                                                           hlp_fn
          write(0,*) '  -- no header found'
          call xerrab ('No header found in spec_pwr.dat file !')
        end if
        backspace(lun)
   30   read(lun,hlp_frm,err=40,end=40) j,(hlp_r,i=1,n_hlp),m
        if(m.eq.1234) go to 30
   40   backspace(lun)
      else

c*** Write the header

        if (ntrgts.eq.2) then
          hlp_s(1) = hlp_c(1)
          hlp_s(2) = hlp_c(4)
        else
          hlp_s = hlp_c
        end if
        rewind(lun)
        write(hlp_frh,'(a,i1,a,i3,a)')
     .   '(a5,14x,',m_hlpa,'a4,7x,',ntrgts*m_hlp,'(a4,a1,7x))'
c        print *,'hlp_frh: ',hlp_frh
        write(lun,hlp_frh) '#case',(hlp_hda(i),i=1,m_hlpa),
     -    ((hlp_hdr(i),trim(hlp_s(j)),i=1,m_hlp),j=1,ntrgts)
      end if

c*** Extract the case number

      call ipgeti('b2mndr_pulse_number', shot )
      if(shot.eq.0) call ipgeti('b2mndr_shot_number', shot )
      call ipgeti('b2mndr_run_number', run )
#if AL_MAJOR_VERSION > 4
      call xertst( 0.le.shot, 'Invalid pulse number')
      call xertst( 0.le.run, 'Invalid run number')
#else
      call xertst( 0.le.shot.and.shot.le.214748, 'Invalid shot number')
      call xertst( 0.le.run.and.run.le.99999, 'Invalid run number')
#endif
      not_default = shot.gt.0.and.run.ge.0
      if (shot.gt.0) then
        write(case_number,'(i6,a1,i5.5)') shot,'/',run
      else
        path = ' '
#if AL_MAJOR_VERSION > 4
        call ipgetc('b2mndr_ids_path', path )
#endif
        call strip_spaces(path)
        not_default = not_default.or..not.streql(path,' ')
        filename='shotnumber.history'
        call find_file(filename,ex)
        if (ex) then
          open(99,file=filename)
   11     continue
          read(99,'(a)',end=12) zeile
          goto 11
        else
          zeile = ' '
        end if
   12   continue
        if (ex) close(99)
        call strip_spaces(zeile)
        if (.not.streql(zeile,' ').or..not.streql(path,' ')) then
          username = usrnam()
          call ipgetc('b2mndr_user', username )
          call xertst( .not.streql(username,' '),
     .     'User name not defined !')
          not_default = not_default.or..not.streql(username,usrnam())
          home_dir = '/home/'//trim(username)
          database = 'solps-iter'
#ifdef IMAS
          write(version,'(i1)') IMAS_MAJOR_VERSION
#else
          version = '0'
#endif
          device_env = ' '
#ifndef NO_GETENV
#ifdef USE_PXFGETENV
          CALL PXFGETENV ('HOME', 0, home_dir, lenval, ierror)
          CALL PXFGETENV ('DEVICE', 0, device_env, lenval, ierror)
#else
          call get_environment_variable('HOME', status=ierror,
     .                                          length=lenval)
          if (ierror.eq.0) call get_environment_variable('HOME',
     .                                          value=home_dir)
          call get_environment_variable('DEVICE', status=ierror,
     .                                          length=lenval)
          if (ierror.eq.0) call get_environment_variable('DEVICE',
     .                                          value=device_env)
#endif
          if (.not.streql(device_env,' ')) database = device_env
          if (streql(database,'iter')) database = 'ITER'
#endif
          call ipgetc('b2mndr_device', database )
          call ipgetc('b2mndr_database', database )
          call xertst( .not.streql(database,' '),
     .     'Database not defined !')
          not_default = not_default.or.
     &                 (.not.streql(database, device_env).and.
     &                  .not.streql(database,'solps-iter'))
          imasdir = trim(home_dir)//'/public/imasdb/'
     &            //trim(database)//'/'//trim(version)
#if AL_MAJOR_VERSION < 5
     &            //'/'//int2str(run/10000)
#endif
#ifndef NO_GETENV
#ifdef USE_PXFGETENV
          CALL PXFGETENV ('IMASDIR', 0, imasdir, lenval, ierror)
#else
          call get_environment_variable('IMASDIR',
     .     status=ierror, length=lenval)
          if (ierror.eq.0)
     >     call get_environment_variable('IMASDIR',value=imasdir)
#endif
#endif
          if (.not.not_default) path = trim(zeile)
          if (streql(path,' ')) path = trim(imasdir)
          l=index(path,'$HOME')
          if (l.gt.0) then
            path = trim(home_dir)//trim(path(l+5:256))
          end if
          l=index(path,'$SOLPSTOP')
          if (l.gt.0) then
            path = trim(solpstop)//trim(path(l+9:256))
          end if
          l=index(path,'$IMASDIR')
          if (l.gt.0) then
            path = trim(imasdir)//trim(path(l+8:256))
          end if
          do while (path(len_trim(path):len_trim(path)).eq.'/')
            path(len_trim(path):len_trim(path))=' '
          end do
          l=index(path,trim(imasdir))
#if AL_MAJOR_VERSION > 4
          if (l.gt.0) then
            m=index(path(l+len_trim(imasdir):256),'/')
            case_number =
     .       trim(path(m+l+len_trim(imasdir):m+l+len_trim(imasdir)+11))
          else
#endif
            l=index(path,'/')
            m=index(path(l+1:len_trim(path)),'/')
            do while (m.ne.0)
              l=l+m
              m=index(path(l+1:len_trim(path)),'/')
            end do
            case_number = path(l+1:l+12)
#if AL_MAJOR_VERSION > 4
          end if
#endif
        end if
        if (streql(case_number,' ')) then
          lv_hlp=0
          l=index(label,'#')
          m=index(label,'_')
          if(m.gt.l+1) then
            hlp_inp=label(l+1:m-1)
            read(hlp_inp,'(i8)',err=99) lv_hlp
          end if
   99     continue
          write(case_number,'(i4)') lv_hlp
        end if
      end if

c*** Write the data

      write(lun,hlp_frm) case_number,(eqa_hlp(j),j=1,m_hlpa),
     ,                 ((equ_hlp(i,j),j=1,m_hlp),i=1,ntrgts),1234
      close(lun)
      return
c======================================================================
  900 write(0,*) 'b2paspwr: unable to open ',hlp_fn
      call xerrab ('Unable to open spec_pwr.dat file!')
c======================================================================
      end subroutine b2paspwr
