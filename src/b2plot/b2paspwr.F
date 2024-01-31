      subroutine b2paspwr(nx,ny,ns,jsep,ic)
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
      use b2mod_types
      use b2mod_geo
      use b2mod_diag
      use b2mod_plasma
      use b2mod_ipmain
      use b2mod_sources
      use b2mod_indirect
      use b2mod_solpstop
      implicit none
c
      integer nx,ny,ns,jsep,ic
      integer lun,m_hlp,m_hlpa
      parameter (lun=96, m_hlp=9, m_hlpa=1)
      integer, save :: shot = 0, run = 0
      character*12 case_number
      character*16 username
      character*24 database,version
      character*256 filename,path,imasdir,home_dir,zeile
      integer i,j,k,l,m,n,ix,ixx,isep,lv_hlp,n_hlp
      real (kind=R8) ::
     ,  hlp_pwr_in(2),hlp_rad_div(2),hlp_rad_trg(2),hlp_area(2),
     ,  hlp_spwr(2),hlp_spwt(2),hlp_rdfr(2),hlp_rad_tot(2),
     ,  hlp_area_t(2)
      real (kind=R8) ::
     ,  equ_hlp(2,m_hlp),hlp_r,s,eqa_hlp(m_hlpa),hlp_tot
      equivalence
     ,  (hlp_tot,eqa_hlp(1))
      equivalence
     ,  (hlp_pwr_in,equ_hlp(1,1)),  (hlp_rad_div,equ_hlp(1,2)),
     ,  (hlp_rad_trg,equ_hlp(1,3)), (hlp_area,equ_hlp(1,4)),
     ,  (hlp_spwr,equ_hlp(1,5)),    (hlp_spwt,equ_hlp(1,6)),
     ,  (hlp_rdfr,equ_hlp(1,7)),    (hlp_rad_tot,equ_hlp(1,8)),
     ,  (hlp_area_t,equ_hlp(1,9))
      logical ex, not_default
      character hlp_fn*24,hlp_c(2),hlp_u,hlp_frm*256,hlp_frh*256,
     ,  hlp_hdr(m_hlp)*4,hlp_hda(m_hlpa)*4,hlp_inp*8
      data hlp_fn / 'spec_pwr.dat' /
      data hlp_c / 'i', 'o' /
      data hlp_hdr / 'Pin_','Pdv_','Ptg_','Sdv_','Psw_',
     ,               'Pst_','RdF_','Prd_','Sdt_' /
      data hlp_hda / 'Pin' /
      integer lnblnk
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

      n_hlp=2*m_hlp+m_hlpa
      isep=jsep+1

c*** Calculate the input power for each divertor

      hlp_tot=0.0_R8
      ix=leftcut(1)-1
      s=-1.0_R8
      do l=1,2
        hlp_pwr_in(l)=0.0_R8
        do j=0,ny-1
          if (l.eq.1) then
            ixx=rightix(ix,j)
          else
            ixx=ix
          end if
          hlp_pwr_in(l)=hlp_pwr_in(l)+fhe(ixx,j,0)+fhi(ixx,j,0)
        end do
        hlp_pwr_in(l)=hlp_pwr_in(l)*s
        hlp_tot=hlp_tot+hlp_pwr_in(l)
        ix=rightcut(1)
        s=1.0_R8
        hlp_pwr_in(l)=1.0e-6_R8*hlp_pwr_in(l)
      end do
      hlp_tot=1.e-6*hlp_tot
c      print *,'Input powers:',hlp_tot,hlp_pwr_in

c*** Calculate the areas and power radiated in the divertors

      m=ic
      n=leftcut(1)-1
      do l=1,2
        hlp_area(l)=0.0_R8
        hlp_area_t(l)=0.0_R8
        hlp_rad_div(l)=0.0_R8
        do i=m,n
          hlp_area(l)=hlp_area(l)+gs(i,isep,1)+gs(i,isep-1,1)
          do j=0,ny-1
            do k=0,ns-1
              hlp_rad_div(l)=hlp_rad_div(l)-rqbrm(i,j,k)
              hlp_rad_div(l)=hlp_rad_div(l)-rqrad(i,j,k)
            end do
          end do
        end do
        hlp_area(l)=0.5_R8*hlp_area(l)
        hlp_rad_div(l)=1.0e-6_R8*hlp_rad_div(l)
        m=rightcut(1)
        n=nx-1-ic
      end do
c      print *,'Areas: ',hlp_area
c      print *,'Div. rad.: ',hlp_rad_div

c*** calculate the power radiated near the targets

      m=0
      n=ic-1
      do l=1,2
        hlp_area_t(l)=hlp_area(l)
        hlp_rad_trg(l)=0.0_R8
        do i=m,n
          hlp_area_t(l)=hlp_area_t(l)+
     .     0.5_R8*(gs(i,isep,1)+gs(i,isep-1,1))
          do j=0,ny-1
            do k=0,ns-1
              hlp_rad_trg(l)=hlp_rad_trg(l)-rqbrm(i,j,k)
              hlp_rad_trg(l)=hlp_rad_trg(l)-rqrad(i,j,k)
            end do
          end do
        end do
        hlp_rad_trg(l)=1.0e-6_R8*hlp_rad_trg(l)
        m=nx-ic
        n=nx-1
      end do
c      print *,'Trg. rad.: ',hlp_rad_trg

c*** total radiated power and its fraction

      do i=1,2
        hlp_rad_tot(i)=hlp_rad_div(i)+hlp_rad_trg(i)
        if(hlp_pwr_in(i).gt.0.0_R8) then
          hlp_rdfr(i)=hlp_rad_tot(i)/hlp_pwr_in(i)
        else
          hlp_rdfr(i)=0.0_R8
        end if
      end do

c*** calculate the specific powers

      do l=1,2
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

        rewind(lun)
        write(hlp_frh,'(a,i1,a,i3,a)')
     .   '(a5,14x,',m_hlpa,'a4,7x,',2*m_hlp,'(a4,a1,7x))'
c        print *,'hlp_frh: ',hlp_frh
        write(lun,hlp_frh) '#case',(hlp_hda(i),i=1,m_hlpa),
     -    ((hlp_hdr(i),hlp_c(j),i=1,m_hlp),j=1,2)
      end if

c*** Extract the case number

      call ipgeti('b2mndr_shot_number', shot )
      call ipgeti('b2mndr_run_number', run )
#if AL_MAJOR_VERSION > 4
      call xertst( 0.le.shot, 'Invalid shot number')
      call xertst( 0.le.run, 'Invalid run number')
#else
      call xertst( 0.le.shot.and.shot.le.214748, 'Invalid shot number')
      call xertst( 0.le.run.and.run.le.99999, 'Invalid run number')
#endif
      if (shot.gt.0) then
        write(case_number,'(i6,a1,i5.5)') shot,'/',run
      else
        path = ' '
#if AL_MAJOR_VERSION > 4
        call ipgetc('b2mndr_ids_path', path )
#endif
        call strip_spaces(path)
        not_default = .not.streql(path,' ')
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
            path = trim(solpstop())//trim(path(l+9:256))
          end if
          l=index(path,'$IMASDIR')
          if (l.gt.0) then
            path = trim(imasdir)//trim(path(l+8:256))
          end if
          do while (path(lnblnk(path):lnblnk(path)).eq.'/')
            path(lnblnk(path):lnblnk(path))=' '
          end do
          l=index(path,trim(imasdir))
#if AL_MAJOR_VERSION > 4
          if (l.gt.0) then
            m=index(path(l+lnblnk(imasdir):256),'/')
            case_number =
     .       trim(path(m+l+lnblnk(imasdir):m+l+lnblnk(imasdir)+11))
          else
#endif
            l=index(path,'/')
            m=index(path(l+1:lnblnk(path)),'/')
            do while (m.ne.0)
              l=l+m
              m=index(path(l+1:lnblnk(path)),'/')
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
     ,                            ((equ_hlp(i,j),j=1,m_hlp),i=1,2),1234
      close(lun)
      return
c======================================================================
  900 write(0,*) 'b2paspwr: unable to open ',hlp_fn
      call xerrab ('Unable to open spec_pwr.dat file!')
c======================================================================
      end subroutine b2paspwr
