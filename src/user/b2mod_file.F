      module b2mod_file
      use b2mod_types
      use b2mod_diag
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none

c======================================================================
c***  Writes the tracing data on the output files for plotting
c***
c***  It can be called by the following entries:
c***
c***    b2filei     to Initialise the files
c***    b2filew     to Write the data on each internal iteration
c***    b2filef     to Flush the buffers
c***
c*** The files are connected to the units
c***
c***    tracing/residuals.trc   to  ntrc(ntrc_res)
c***    tracing/sources.trc     to  ntrc(ntrc_src)
c----------------------------------------------------------------------
c
c*** lvr, lvs:    version numbers to appear in the files
c
      integer lvr, lvs, lvt
c      parameter (lvr=960719, lvs=960821, lvt=980805)
      parameter (lvr=960719, lvs=960821, lvt=20001218)
c
c*** nf_tmp: number of files which are written in b2filew
c*** lf_tmp: number of files whose buffers should be flushed in b2filef
c***        (must be consistent with ntrc assignment)
c***
c***    Current assignment is (see b2mod_diag.F):
c***        ntrc(ntrc_res)     'tracing/residuals.trc'
c***        ntrc(ntrc_src)     'tracing/sources.trc'
c***        ntrc(ntrc_tst)     'tracing/test.trc'
c
      integer :: nh_tmp, nf_tmp
      parameter (nh_tmp=6, nf_tmp=3)
      integer, save :: nis_tmp
c
      character, save :: hlp_frmr*160, hlp_frms*160, hlp_frmt*160,
     ,          hlp_emsg(5)*30
      data hlp_emsg /
     1      'cannot open the file',
     2      'wrong format of the file',
     3      'unexpected error. Unopened',
     4      'error during writing -->',
     5      'unexpected error. Cannot close' /

      contains
c
      subroutine b2filei(ns,nregionv)
      use b2mod_plasma
      use b2mod_version
      implicit none
c
      integer ns, nregionv
c
      integer, save :: i, j, k, l, m, n, lf_tmp
      integer :: mh_tmp(0:ns-1)
      logical, save :: ex
      character, save :: hlp_c, f*(lenfnm), hlp_s(3)
      character :: hlp_h(-1:nh_tmp)*2, hlp_al(0:ns-1)*2
      real (kind=R8), save :: stim = 0.0_R8
c
c----------------------------------------------------------------------
      data hlp_h / 'ms', 'es', 'ee', 'ei', 'co', 'mo', 'ni', 'mo' /
      data hlp_s / 'r', 'c', 'i' /
c======================================================================
c
c----------------------------------------------------------------------
c*** Initialisation call.
c*** Check whether the directory "tracing" exists and the files can
c*** be opened. If yes, then the data will be written to the files.
c*** The file headers are written to the files if the files are empty.
c----------------------------------------------------------------------
      call subini('b2filei')
c
      write(*,*) 'b2file: nstrai, nnatmi = ', nstrai, nnatmi
      lf_tmp=ntrc_used
c
c*** check the number of units specified in ntrc
c
      call xertst(lf_tmp.le.len_ntrc,'b2file: too many tracing units')
      call ipgetr('b2mndr_stim', stim)
c----------------------------------------------------------------------
c*** Test data - temporary
c
      if(.not.filedata(ntrc_tst)) go to 1000
c
      filedata(ntrc_tst)=.false.
      f=fn_ntrc(ntrc_tst)
      call find_file(f,ex)
      l=1
      open(ntrc(ntrc_tst),file=f,err=900)
      rewind(ntrc(ntrc_tst))
      if(ex.and.stim.lt.0.0_R8) then
c*** Old file to be appended - check it
        l=2
        read(ntrc(ntrc_tst),*,err=900) hlp_c,i,m,n
        if(i.ne.lvt) then
          write(*,'(a,i8,3a,i8)')
     ,     '!! b2file: current version ',lvt,
     ,     ' differs from that in the file ',trim(f),' : ',i
          write(*,*) ' --  Remove the file first!'
          go to 1000
        end if
        if(m.ne.nstrai) then
          write(*,'(a,i3,3a,i3)')
     ,     '!! b2file: current number of strata (nstrai=',nstrai,
     ,     ') differs from that in the file ',trim(f),' : ',m
          write(*,*) ' --  Remove the file first!'
          go to 1000
        end if
        if(n.ne.nnatmi) then
          write(*,'(a,i2,3a,i2)')
     ,     '!! b2file: current number of atoms (nnatmi=',nnatmi,
     ,     ') differs from that in the file ',trim(f),' : ',n
          write(*,*) ' --  Remove the file first!'
          go to 1000
        end if
      else
c*** New file - write the headers
        l=4
        write(ntrc(ntrc_tst),700,err=900) '#',lvt,nstrai,nnatmi
        write(ntrc(ntrc_tst),710,err=900) '## ',trim(label)

        if(nstrai.gt.0) then
          if (nnatmi.gt.0) then
            write(hlp_frmt,'(a,i3,a,i2,a)')
     /         '(''# it'',7x,''time'',6x,',nstrai*5,
     /         '(1x,a8,i2.2,a3,',nnatmi,'(a8,i2.2,a1,a2,1x)))'
          else
            write(hlp_frmt,'(a,i3,a,i2,a)')
     /         '(''# it'',7x,''time'',6x,',nstrai*5,
     /         '(1x,a8,i2.2,a3))'
          end if
#ifdef DBG
          print *,'%%% b2file: header format'   !###
          print *,trim(hlp_frmt)                !###
#endif
          if (nnatmi.gt.0) then
            write(ntrc(ntrc_tst),hlp_frmt,err=900)
     -       ('sclstra_',i,'   ',
     ,       ('sclstra_',i,'_',atoms_list(k),k=1,nnatmi),
     ,        'sclalph_',i,'   ',
     ,       ('sclalph_',i,'_',atoms_list(k),k=1,nnatmi),
     ,        'sclbase_',i,'   ',
     ,       ('sclbase_',i,'_',atoms_list(k),k=1,nnatmi),
     ,        'sclrtio_',i,'   ',
     ,       ('sclrtio_',i,'_',atoms_list(k),k=1,nnatmi),
     ,        'sclsrcn_',i,'   ',
     ,       ('sclsrcn_',i,'_',atoms_list(k),k=1,nnatmi),i=1,nstrai)
          else
            write(ntrc(ntrc_tst),hlp_frmt,err=900)
     -       ('sclstra_',i,'   ',
     ,        'sclalph_',i,'   ',
     ,        'sclbase_',i,'   ',
     ,        'sclrtio_',i,'   ',
     ,        'sclsrcn_',i,'   ',i=1,nstrai)
          end if
        end if
      end if
      if (nstrai.gt.0) write(hlp_frmt,'(a,i3,a)')
     -      '(1x,i3.3,1p,e16.7,',nstrai*5*(nnatmi+1),'(e12.4,2x))'
      filedata(ntrc_tst)=.true.
      go to 1000
c
  900 write(*,'(1x,2a,1x,a)') 'b2file: ',trim(hlp_emsg(l)),trim(f)
c----------------------------------------------------------------------
c*** Residuals
c
 1000 continue
      if(.not.filedata(ntrc_res)) go to 10
c
      filedata(ntrc_res)=.false.
      f=fn_ntrc(ntrc_res)
      call find_file(f,ex)
      l=1
      open(ntrc(ntrc_res),file=f,err=5)
      rewind(ntrc(ntrc_res))
      if(ex.and.stim.lt.0.0_R8) then
        l=2
        read(ntrc(ntrc_res),*,err=5) hlp_c,i,j,k
        if(i.ne.lvr) then
          write(*,'(a,i8,3a,i8)')
     ,     '!! b2file: current version ',lvr,
     ,     ' differs from that in the file ',trim(f),' : ',i
          write(*,*) ' --  Remove the file first!'
          go to 10
        end if
        if(j.ne.nregionv) then
          write(*,'(a,i2,3a,i2)')
     ,     '!! b2file: current number of regions ',nregionv,
     ,     ' differs from that in the file ',trim(f),' : ',j
          write(*,*) ' --  Remove the file first!'
          go to 10
        end if
        if(k.ne.ns-nspecies) then
          write(*,'(a,i3,3a,i3)')
     ,     '!! b2file: current number of ion species (nis=',ns-nspecies,
     ,     ') differs from that in the file ',trim(f),' : ',k
          write(*,*) ' --  Remove the file first!'
          go to 10
        else
          nis_tmp=k
        end if
      else
c*** New file - write the headers
        k=0
        do j=1,nspecies
          do i=1,nfluids(j)
            k=k+1
            mh_tmp(k)=i
            hlp_al(k)=species_list(j)
          end do
        end do
        l=4
        nis_tmp=k

        write(ntrc(ntrc_res),700,err=5) '#',lvr,nregionv,nis_tmp
        write(ntrc(ntrc_res),710,err=5) '## ',trim(label)

        write(hlp_frmr,'(a,i3,a,i3,a)')
     /        '(''# it'',7x,''time'',8x,',
     -         nregionv*2,
     /        '(2x,a3,a2,i2.2,5x),',
     -         nregionv*(nis_tmp)*2,
     /        '(a3,a2,i2.2,''_'',i2.2,''_'',a2,1x))'
#ifdef DBG
        print *,'%%% b2file: header format'   !###
        print *,trim(hlp_frmr)                !###
#endif
        write(ntrc(ntrc_res),hlp_frmr,err=5)
     ,    (('res',hlp_h(i),j,j=1,nregionv),i=1,2),
     ,   ((('res',hlp_h(i),j,mh_tmp(k),hlp_al(k),j=1,nregionv),
     ,                                           i=3,4),k=1,nis_tmp)
      end if
      write(hlp_frmr,'(a,i3,a)')
     - '(1x,i3.3,1p,',nregionv*2*(nis_tmp+1),'(2x,e12.4))'
      filedata(ntrc_res)=.true.
      go to 10
c
    5 write(*,'(1x,2a,1x,a)') 'b2file: ',trim(hlp_emsg(l)),trim(f)
c----------------------------------------------------------------------
c*** Sources
c
   10 if(.not.filedata(ntrc_src)) go to 35
c
      filedata(ntrc_src)=.false.
      f=fn_ntrc(ntrc_src)
      call find_file(f,ex)
      l=1
      open(ntrc(ntrc_src),file=f,err=30)
      rewind(ntrc(ntrc_src))
      if(ex.and.stim.lt.0.0_R8) then
        l=2
        read(ntrc(ntrc_src),700,err=30) hlp_c,i,j,k
        if(i.ne.lvs) then
          write(*,'(a,i8,3a,i8)')
     ,     '!! b2file: current version ',lvs,
     ,     ' differs from that in the file ',trim(f),' : ',i
          write(*,*) ' --  Remove the file first!'
          go to 35
        end if
        if(j.ne.nregionv) then
          write(*,'(a,i2,3a,i2)')
     ,     '!! b2file: current number of regions ',nregionv,
     ,     ' differs from that in the file ',trim(f),' : ',j
          write(*,*) ' --  Remove the file first!'
          go to 35
        end if
        if(k.ne.nspecies) then
          write(*,'(a,i2,3a,i2)')
     ,     '!! b2file: current number of species (nspecies=',nspecies,
     ,     ') differs from that in the file ',trim(f),' : ',k
          write(*,*) ' --  Remove the file first!'
          go to 35
        end if
      else
c*** New file - write the headers
        l=4
        write(ntrc(ntrc_src),700,err=30) '#',lvs,nregionv,nspecies
        write(ntrc(ntrc_src),710,err=30) '## ',trim(label)

        write(hlp_frms,'(a,i3,a,i3,a)')
     /        '(''# it'',2x,',
     -         (nregionv+1)*12,
     /        '(3x,a3,a2,i2.2,a1,3x),',
     -         (nregionv+1)*nspecies*6,
     /        '(2x,a3,a2,i2.2,a1,''_'',a2,1x))'
#ifdef DBG
        print *,'%%% b2file: header format'
        print *,trim(hlp_frms)
#endif
        write(ntrc(ntrc_src),hlp_frms,err=30)
     ,    ((('src',hlp_h(i),j,hlp_s(l),l=1,3),j=0,nregionv),i=-1,2),
     ,   (((('src',hlp_h(i),j,hlp_s(l),species_list(k),l=1,3),
     ,                           j=0,nregionv),i=5,6),k=1,nspecies)
      end if
      write(hlp_frms,'(a,i3,a)')
     - '(1x,i3.3,1p,',(nregionv+1)*6*(nspecies+2),'(2x,e12.4))'
      filedata(ntrc_src)=.true.
      go to 35

   30 write(*,'(1x,2a,1x,a)') 'b2file: ',trim(hlp_emsg(l)),trim(f)
c----------------------------------------------------------------------
c*** Position the files for appending the data
c
   35 do k=1,nf_tmp
        l=k+1
        if(filedata(l)) then
          f=fn_ntrc(l)
          rewind(ntrc(l))
   20     read(ntrc(l),'(a1)',end=25) hlp_c
          go to 20
   25     backspace ntrc(l)
        end if
      end do
      call subend()
      return
c======================================================================
  700 format(a1,2x,4i10)
  710 format(2a)
      end subroutine b2filei
c======================================================================
      subroutine b2filec
      implicit none
      integer :: k, l, iline
      character f*(lenfnm), hlp_c
      logical :: ex
      do k=1,nf_tmp
        l=k+1
        if(filedata(l)) then
          f=fn_ntrc(l)
          call find_file(f,ex)
          open(ntrc(l),file=f,err=290)
          rewind(ntrc(l))
          iline = 0
  200     read(ntrc(l),'(a1)',end=205) hlp_c
          if (hlp_c.eq.'#') then
            iline = iline + 1
          end if
          go to 200
  205     backspace ntrc(l)
        end if
        goto 299
  290   write(*,*) 'b2filec: cannot open file ',trim(f)
        filedata(l)=.false.
        iline=-1
  299   continue
      end do
      return
      end subroutine b2filec
c======================================================================
      subroutine b2filew(it,nregionv)
      implicit none
      integer, intent(in) :: it, nregionv
      integer i, j, k, l
      character f*(lenfnm)
      logical :: op
c
c----------------------------------------------------------------------
c*** Write the data onto the files
c----------------------------------------------------------------------

      if(filedata(ntrc_tst)) then
        f=fn_ntrc(ntrc_tst)
        inquire(ntrc(ntrc_tst),opened=op)
        l=3
        if(.not.op) go to 62
        l=4
        if (nstrai.gt.0) write(ntrc(ntrc_tst),hlp_frmt,err=62) it, tim,
     ,    ((sclstra(i,k),i=0,nnatmi),(sclalph(i,k),i=0,nnatmi),
     ,     (sclbase(i,k),i=0,nnatmi),(sclrtio(i,k),i=0,nnatmi),
     ,     (sclsrcn(i,k),i=0,nnatmi),k=1,nstrai)
        go to 65
   62   write(*,'(1x,2a,1x,a)') 'b2file: ',trim(hlp_emsg(l)),trim(f)
        filedata(ntrc_tst)=.false.
      end if
c
   65 if(filedata(ntrc_res)) then
        f=fn_ntrc(ntrc_res)
        inquire(ntrc(ntrc_res),opened=op)
        l=3
        if(.not.op) go to 72
        l=4
        write(ntrc(ntrc_res),hlp_frmr,err=72) it,
     ,     (reshereg(i),i=1,nregionv),
     ,     (reshireg(i),i=1,nregionv),
     ,    ((rescoreg(i,ions_list(j)),i=1,nregionv),
     ,     (resmoreg(i,ions_list(j)),i=1,nregionv),j=1,nis_tmp)
        go to 75
   72   write(*,'(1x,2a,1x,a)') 'b2file: ',trim(hlp_emsg(l)),trim(f)
        filedata(ntrc_res)=.false.
      end if
c
   75 if(filedata(ntrc_src)) then
        f=fn_ntrc(ntrc_src)
        inquire(ntrc(ntrc_src),opened=op)
        l=3
        if(.not.op) go to 77
        l=4
        write(ntrc(ntrc_src),hlp_frms,err=77) it,
     ,          ((src_ms(i,k),k=1,3),i=0,nregionv),
     ,          ((src_es(i,k),k=1,3),i=0,nregionv),
     ,          ((src_ee(i,k),k=1,3),i=0,nregionv),
     ,          ((src_ei(i,k),k=1,3),i=0,nregionv),
     ,         (((src_ni(i,k,j),k=1,3),i=0,nregionv),
     ,          ((src_mo(i,k,j),k=1,3),i=0,nregionv),j=1,nspecies)
        return
   77   write(*,'(1x,2a,1x,a)') 'b2file: ',trim(hlp_emsg(l)),trim(f)
        filedata(ntrc_src)=.false.
      end if
      return
      end subroutine b2filew
c======================================================================
c
      subroutine b2filef
      implicit none
      integer j, k, l
      character f*(lenfnm)
c----------------------------------------------------------------------
c*** flush the buffers
c----------------------------------------------------------------------
c
c
      do 100 k=1,nf_tmp
        j=k+1
        if(filedata(j)) then
          f=fn_ntrc(j)
#ifdef IBMaix
          call flush_(ntrc(j))
#else
          l=5
          close(ntrc(j),err=86)
#ifndef F77
          l=1
          open(ntrc(j),file=f,err=86,position='append')
#else
          l=1
          open(ntrc(j),file=f,err=86)
   80     read(ntrc(j),'(a1)',end=82) hlp_c
          go to 80
   82     backspace ntrc(j)
#endif
#endif
          go to 100
c
   86     write(*,'(1x,2a,1x,a)') 'b2file: ',trim(hlp_emsg(l)),trim(f)
          filedata(j)=.false.
        end if
  100 continue
      return
      end subroutine b2filef
c
      end module b2mod_file
