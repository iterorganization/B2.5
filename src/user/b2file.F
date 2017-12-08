      subroutine b2file(it,ns,nnreg)
      use b2mod_types
      use b2mod_plasma
      use b2mod_version
      use b2mod_b2cmpa
      use b2mod_diag
c
c  version : 21.12.2000 19:05
c
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
      implicit none
c
      integer ns, nnreg(0:2), it
c
      integer, save ::  i, j, k, l, m, n, lf_tmp, nis_tmp
      integer :: mh_tmp(0:ns-1), nh_tmp, nf_tmp, lvr, lvs, lvt
      parameter (nh_tmp=6, nf_tmp=3)
      logical, save :: ex, op
c
c*** lvr, lvs:    version numbers to appear in the files
c
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
      character, save :: tab, hlp_c, f*(lenfnm), hlp_s(3),
     ,          hlp_frmr*160, hlp_frms*160, hlp_frmt*160,
     ,          hlp_emsg(5)*30
      character :: hlp_h(-1:nh_tmp)*2, hlp_al(0:ns-1)*2
      real (kind=R8), save :: stim = 0.0_R8
      integer iline
c
      external subini, subend, ipgetr, xertst
c
c
c----------------------------------------------------------------------
      data hlp_h / 'ms', 'es', 'ee', 'ei', 'co', 'mo', 'ni', 'mo' /
      data hlp_s / 'r', 'c', 'i' /
      data hlp_emsg /
     1      'cannot open the file',
     2      'wrong format of the file',
     3      'unexpected error. Unopened',
     4      'error during writing -->',
     5      'unexpected error. Cannot close' /
c======================================================================
c
      entry b2filei(it,ns,nnreg)
c----------------------------------------------------------------------
c*** Initialisation call.
c*** Check whether the directory "tracing" exists and the files can
c*** be opened. If yes, then the data will be written to the files.
c*** The file headers are written to the files if the files are empty.
c----------------------------------------------------------------------
      call subini('b2filei')
c
      write(*,*) 'b2file: nstrai, nnatmi = ', nstrai, nnatmi
      tab=char(9)
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
        read(ntrc(ntrc_tst),700,err=900) hlp_c,i,m,n,j
        if(i.ne.lvt) then
          write(*,*) '!! b2file: current version ',lvt,
     ,               ' differs from that in the file ',trim(f)
          write(*,*) ' --  Remove the file first!'
          go to 1000
        end if
        if(m.ne.nstrai) then
          write(*,*) '!! b2file: current number of strata (nstrai=',
     ,     nstrai,') differs from that in the file ',trim(f)
          write(*,*) ' --  Remove the file first!'
          go to 1000
        end if
        if(n.ne.nnatmi) then
          write(*,*) '!! b2file: current number of atoms (nnatmi=',
     ,     nnatmi,') differs from that in the file ',trim(f)
          write(*,*) ' --  Remove the file first!'
          go to 1000
        end if
      else
c*** New file - write the headers
        l=4
        write(ntrc(ntrc_tst),700,err=900) '#',lvt,nstrai,nnatmi
        write(ntrc(ntrc_tst),710,err=900) '## ',label

        if(nstrai.gt.0) then
c          write(hlp_frmt,'(7h(3h#it,,i3,21h(1h'//tab//
c     /     ',a8,i2.2,1h_,a2)))')
c     -                          nstrai*5*(nnatmi+1)
          write(hlp_frmt,'(a,i3,a)')
     /         '(''#it'',',
     -                          nstrai*5*(nnatmi+1),
     /         '('''//tab//''',a8,i2.2,''_'',a2))'
c<<<
          print *,'%%% b2file: header format'   !###
          print *,trim(hlp_frmt)                !###
c>>>
          write(ntrc(ntrc_tst),hlp_frmt,err=900)
     -     ('sclstra_',i,'  ',
     ,     ('sclstra_',i,atoms_list(k),k=1,nnatmi),
     ,     'sclalph_',i,'  ',
     ,     ('sclalph_',i,atoms_list(k),k=1,nnatmi),
     ,      'sclbase_',i,'  ',
     ,     ('sclbase_',i,atoms_list(k),k=1,nnatmi),
     ,      'sclrtio_',i,'  ',
     ,     ('sclrtio_',i,atoms_list(k),k=1,nnatmi),
     ,      'sclsrcn_',i,'  ',
     ,     ('sclsrcn_',i,atoms_list(k),k=1,nnatmi),i=1,nstrai)
        end if
      end if
      if (nstrai.gt.0)
     > write(hlp_frmt,'(12h(1x,i3.3,1p,,i3,12h(1h'//tab//',e12.4)))')
     -                          nstrai*5*(nnatmi+1)
      filedata(ntrc_tst)=.true.
      go to 1000
c
 900  write(*,'(2a,1x,a)') 'b2file: ',trim(hlp_emsg(l)),trim(f)
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
        read(ntrc(ntrc_res),700,err=5) hlp_c,i,j,k
        if(i.ne.lvr) then
          write(*,*) '!! b2file: current version ',lvr,
     ,               ' differs from that in the file ',trim(f)
          write(*,*) ' --  Remove the file first!'
          go to 10
        end if
        if(j.ne.nnreg(0)) then
          write(*,*) '!! b2file: current number of regions ',
     ,                nnreg(0),' differs from that in the file ',trim(f)
          write(*,*) ' --  Remove the file first!'
          go to 10
        end if
        if(k.ne.ns-nspecies) then
          write(*,*) '!! b2file: current number of ion species (nis=',
     ,     ns-nspecies,') differs from that in the file ',trim(f)
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

        write(ntrc(ntrc_res),700,err=5) '#',lvr,nnreg(0),nis_tmp
        write(ntrc(ntrc_res),710,err=5) '## ',label

c       write(hlp_frmr,'(7h(3h#it,,i3,17h(1h'//tab//
c     /  ',a3,a2,i2.2),,i3,29h(1h'//tab//
c     /  ',a3,a2,i2.2,1h_,i2.2,a2)))') nnreg(0)*2,nnreg(0)*ns*2
        write(hlp_frmr,'(a,i3,a,i3,a)')
     /        '(''#it'',',
     -         nnreg(0)*2,
     /        '('''//tab//''',a3,a2,i2.2),',
     -         nnreg(0)*(nis_tmp)*2,
     /        '('''//tab//''',a3,a2,i2.2,''_'',i2.2,a2))'
c<<<
        print *,'%%% b2file: header format'   !###
        print *,trim(hlp_frmr)                !###
c>>>

        write(ntrc(ntrc_res),hlp_frmr,err=5)
     ,   (('res',hlp_h(i),j,j=1,nnreg(0)),i=1,2),
     ,   ((('res',hlp_h(i),j,mh_tmp(k),hlp_al(k),j=1,nnreg(0)),
     ,   i=3,4),k=1,nis_tmp)
      end if
      write(hlp_frmr,'(12h(1x,i3.3,1p,,i3,12h(1h'//tab//',e12.4)))')
     -  nnreg(0)*2*(nis_tmp+1)
      filedata(ntrc_res)=.true.
      go to 10
c
  5   write(*,'(2a,1x,a)') 'b2file: ',trim(hlp_emsg(l)),trim(f)
c----------------------------------------------------------------------
c*** Sources
c
 10   if(.not.filedata(ntrc_src)) go to 35
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
          write(*,*) '!! b2file: current version ',lvs,
     ,     ' differs from that in the file ',trim(f)
          write(*,*) ' --  Remove the file first!'
          go to 35
        end if
        if(j.ne.nnreg(0)) then
          write(*,*) '!! b2file: current number of regions ',
     ,                nnreg(0),' differs from that in the file ',trim(f)
          write(*,*) ' --  Remove the file first!'
          go to 35
        end if
        if(k.ne.nspecies) then
          write(*,*) '!! b2file: current number of species (nspecies=',
     ,     nspecies,') differs from that in the file ',trim(f)
          write(*,*) ' --  Remove the file first!'
          go to 35
        end if
      else
c*** New file - write the headers
        l=4
        write(ntrc(ntrc_src),700,err=30) '#',lvs,nnreg(0),nspecies
        write(ntrc(ntrc_src),710,err=30) '## ',label

c        write(hlp_frms,'(7h(3h#it,,i3,20h(1h'//tab//
c     /          ',a3,a2,i2.2,a1),,i3,27h(1h'//tab//
c     /          ',a3,a2,i2.2,a1,1h_,a2)))')
c     -                          (nnreg(0)+1)*12,(nnreg(0)+1)*nspecies*6
        write(hlp_frms,'(a,i3,a,i3,a)')
     /        '(''#it'',',
     -                          (nnreg(0)+1)*12,
     /        '('''//tab//''',a3,a2,i2.2,a1),',
     -                          (nnreg(0)+1)*nspecies*6,
     /        '('''//tab//''',a3,a2,i2.2,a1,''_'',a2))'

c<<<
        print *,'%%% b2file: header format'
        print *,trim(hlp_frms)
c>>>
        write(ntrc(ntrc_src),hlp_frms,err=30)
     ,   ((('src',hlp_h(i),j,hlp_s(l),l=1,3),j=0,nnreg(0)),i=-1,2),
     ,   (((('src',hlp_h(i),j,hlp_s(l),species_list(k),l=1,3),
     ,   j=0,nnreg(0)),i=5,6),k=1,nspecies)
      end if
      write(hlp_frms,'(12h(1x,i3.3,1p,,i3,12h(1h'//tab//',e12.4)))')
     -                                      (nnreg(0)+1)*6*(nspecies+2)

      filedata(ntrc_src)=.true.
      go to 35

 30   write(*,'(2a,1x,a)') 'b2file: ',trim(hlp_emsg(l)),trim(f)
c----------------------------------------------------------------------
c*** Position the files for appending the data
c
 35   do k=1,nf_tmp
        l=k+1
        if(filedata(l)) then
          f=fn_ntrc(l)
          rewind(ntrc(l))
 20       read(ntrc(l),'(a1)',end=25) hlp_c
          go to 20
 25       backspace ntrc(l)
        end if
      end do
      call subend()
      return
c======================================================================
      entry b2filec(it,ns,nnreg)
      do k=1,nf_tmp
        l=k+1
        if(filedata(l)) then
          f=fn_ntrc(l)
          call find_file(f,ex)
          open(ntrc(l),file=f,err=290)
          rewind(ntrc(l))
          iline = 0
 200      read(ntrc(l),'(a1)',end=205) hlp_c
          if (hlp_c.eq.'#') then
            iline = iline + 1
          end if
          go to 200
 205      backspace ntrc(l)
        end if
        goto 299
 290    write(*,*) 'b2filec: cannot open file ',trim(f)
        filedata(l)=.false.
        iline=-1
 299    continue
      end do
      return
c======================================================================
      entry b2filew(it,ns,nnreg)
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
        if (nstrai.gt.0) write(ntrc(ntrc_tst),hlp_frmt,err=62) it,
     ,    ((sclstra(i,k),i=0,nnatmi),(sclalph(i,k),i=0,nnatmi),
     ,     (sclbase(i,k),i=0,nnatmi),(sclrtio(i,k),i=0,nnatmi),
     ,     (sclsrcn(i,k),i=0,nnatmi),k=1,nstrai)
        go to 65
 62     write(*,'(2a,1x,a)') 'b2file: ',trim(hlp_emsg(l)),trim(f)
        filedata(ntrc_tst)=.false.
      end if
c
 65   if(filedata(ntrc_res)) then
        f=fn_ntrc(ntrc_res)
        inquire(ntrc(ntrc_res),opened=op)
        l=3
        if(.not.op) go to 72
        l=4
        write(ntrc(ntrc_res),hlp_frmr,err=72) it,
     ,     (reshereg(i),i=1,nnreg(0)),
     ,     (reshireg(i),i=1,nnreg(0)),
     ,    ((rescoreg(i,ions_list(j)),i=1,nnreg(0)),
     ,     (resmoreg(i,ions_list(j)),i=1,nnreg(0)),j=1,nis_tmp)
        go to 75
 72     write(*,'(2a,1x,a)') 'b2file: ',trim(hlp_emsg(l)),trim(f)
        filedata(ntrc_res)=.false.
      end if
c
 75   if(filedata(ntrc_src)) then
        f=fn_ntrc(ntrc_src)
        inquire(ntrc(ntrc_src),opened=op)
        l=3
        if(.not.op) go to 77
        l=4
        write(ntrc(ntrc_src),hlp_frms,err=77) it,
     ,          ((src_ms(i,k),k=1,3),i=0,nnreg(0)),
     ,          ((src_es(i,k),k=1,3),i=0,nnreg(0)),
     ,          ((src_ee(i,k),k=1,3),i=0,nnreg(0)),
     ,          ((src_ei(i,k),k=1,3),i=0,nnreg(0)),
     ,         (((src_ni(i,k,j),k=1,3),i=0,nnreg(0)),
     ,          ((src_mo(i,k,j),k=1,3),i=0,nnreg(0)),j=1,nspecies)
        return
 77     write(*,'(2a,1x,a)') 'b2file: ',trim(hlp_emsg(l)),trim(f)
        filedata(ntrc_src)=.false.
      end if
      return
c======================================================================
c
      entry b2filef(it,ns,nnreg)
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
 80       read(ntrc(j),'(a1)',end=82) hlp_c
          go to 80
 82       backspace ntrc(j)
#endif
#endif
          go to 100
c
 86       write(*,'(2a,1x,a)') 'b2file: ',trim(hlp_emsg(l)),trim(f)
          filedata(j)=.false.
        end if
 100  continue
      return
c======================================================================
  700 format(a1,2x,4i10)
  710 format(2a)
c
      end subroutine b2file
