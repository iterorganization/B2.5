      subroutine b2trace(nput,nx,ny,ns,nxtl,nxtr,ismain,dtim,BoRiS)
      use b2mod_types
      use b2mod_diag
      use b2mod_b2cmpa
      use b2mod_indirect
      use b2mod_user_namelist
      use b2mod_neutrals_namelist
#ifdef B25_EIRENE
      use eirmod_extrab25
#endif
c
c  version : 29.08.2000 00:26
c
c======================================================================
c*** Output of tracing information at the end of each outer iteration
c*** If an error occurs while opening the files, the balance output is
c*** redirected to the log file for the first call and there is no
c*** further output
c----------------------------------------------------------------------
      implicit none
      integer nput, nx, ny, ns, nxtl, nxtr, ismain
      real (kind=R8) :: dtim, BoRiS
      integer iline
      integer, save :: use_recyceir
      logical, save :: tracing_initialised
      data tracing_initialised /.false./
      data use_recyceir /1/
#ifdef B25_EIRENE
      integer ix, iy, is, istra
      real (kind=R8) :: recyctmp
c----------------------------------------------------------------------
c
c*** Correct the neutral fluxes according to recyc values
c
      write(*,'(a,6i5)') 'nspecies,nnatmi,nnstsi,nnlimi,nnmoli,nnplsi=',
     1 nspecies,nnatmi,nnstsi,nnlimi,nnmoli,nnplsi
csw 2012 nlmpgs --> nnstsi+nnlimi (known after 1st EIRENE call)
      do ix=1,nnstsi+nnlimi
        wldnek(ix,0)=0.0_R8
        wldnep(ix,0)=0.0_R8
        do is=1,nnatmi
          wldna(ix,is,0)=0.0_R8
          wldra(ix,is,0)=0.0_R8
          ewlda(ix,is,0)=0.0_R8
        end do
        do is=1,nnmoli
          wldnm(ix,is,0)=0.0_R8
          wldrm(ix,is,0)=0.0_R8
          ewldm(ix,is,0)=0.0_R8
        end do
        if((ix.le.nnlimi).or.
     +     (ix.gt.nnlimi.and.ix.le.nnlimi+nnstsi)) then
          do istra=1,nstrai
            if(use_recyceir.gt.0) then
              recyctmp=recyceir(istra)
            else
              recyctmp=recyc(ismain,istra)
            endif
            wldnek(ix,0)=wldnek(ix,0)+
     +       wldnek(ix,istra)*recyctmp
            wldnep(ix,0)=wldnep(ix,0)+
     +       wldnep(ix,istra)*recyctmp
            do is=1,nnatmi
              iy=eb2atcr(is)
              if(use_recyceir.eq.0) recyctmp=recyc(iy,istra)
              wldna(ix,is,0)=wldna(ix,is,0)+
     +         wldna(ix,is,istra)*recyctmp
              wldra(ix,is,0)=wldra(ix,is,0)+
     +         wldra(ix,is,istra)*recyctmp
              ewlda(ix,is,0)=ewlda(ix,is,0)+
     +         ewlda(ix,is,istra)*recyctmp
            end do
            do is=1,nnmoli
              if (lmolscl(is).gt.0 .and. lmolscl(is).le.nnatmi) then
                iy=eb2atcr(lmolscl(is))
                if(use_recyceir.eq.0) recyctmp=recyc(iy,istra)
                wldnm(ix,is,0)=wldnm(ix,is,0)+
     +           wldnm(ix,is,istra)*recyctmp
                wldrm(ix,is,0)=wldrm(ix,is,0)+
     +           wldrm(ix,is,istra)*recyctmp
                ewldm(ix,is,0)=ewldm(ix,is,0)+
     +           ewldm(ix,is,istra)*recyctmp
              else if (use_recyceir.gt.0) then !here we are rescaling by stratum, not by species
                wldnm(ix,is,0)=wldnm(ix,is,0)+
     +           wldnm(ix,is,istra)*recyctmp
                wldrm(ix,is,0)=wldrm(ix,is,0)+
     +           wldrm(ix,is,istra)*recyctmp
                ewldm(ix,is,0)=ewldm(ix,is,0)+
     +           ewldm(ix,is,istra)*recyctmp
              else
                if (ix.eq.1 .and. istra.eq.1) then
                  write(*,*) 'Faulty or missing LMOLSCL assignment:'
                  write(*,*) 'LMOLSCL(',is,')=',lmolscl(is)
                  write(*,*) 'Should be valued between 1 and NATM'
                  write(*,*) 'Modify in b2.neutrals.parameters'
                  write(*,*) 'See manual or b2cdcn.f file for usage'
                  write(*,*) 'Scaling of molecular tallies skipped'
                endif
                wldnm(ix,is,0)=wldnm(ix,is,0)+wldnm(ix,is,istra)
                wldrm(ix,is,0)=wldrm(ix,is,0)+wldrm(ix,is,istra)
                ewldm(ix,is,0)=ewldm(ix,is,0)+ewldm(ix,is,istra)
              endif
            end do
          end do
        endif
      end do
      if (nnatmi.ne.nspecies) write (*,"(a,2i4)")
     &  'b2trace: partial source rescaling: natmi,nspecies=',
     &                                     nnatmi,nspecies  !srv 18.05.09
#endif
      if(ank_tracing.eq.0) return
c*** Particle balance in old style by Andrei Kukushkin
      if(filedata(ntrc_bln)) call b2blnc(nx,ny,ns,dtim,iline)
c*** Particle balance in SPb style (for use in the Method of effective source for code speed-up,
c***  see E. Kaveeva et al 2018 Nucl. Fusion 58 126018)
      if(filedata(ntrc_bls)) call b2blns(nx,ny,ns,dtim,iline)
c*** Power balance
      if(filedata(ntrc_ble)) call b2blne(nx,ny,ns,BoRiS,iline)
c*** Parallel momentum balance
      if(filedata(ntrc_blm)) call b2blnm(nx,ny,ns,iline)
c*** Additional integrals
      if(filedata(ntrc_int)) call b2wrint(nx,ny,ns,nnreg(0),iline)
      if(filedata(ntrc_srt)) call b2wrsrt(nx,ny,ns,nnreg(0),iline)
c*** Separatrix data
      if(filedata(ntrc_sep)) call b2wrsep(nx,ny,ns,BoRiS,iline)
c*** User-specific data
      if(filedata(ntrc_usr))
     > call b2usrtrc(nx,ny,ns,nxtl,nxtr,BoRiS,iline)
c
      return
c======================================================================
      entry b2trci(nput,nx,ny,ns,nxtl,nxtr,ismain,dtim,BoRiS)
c*** Initialisation call
      if (tracing_initialised) return
      call alloc_b2mod_user(nx,ns,nlim,nmol,ntns)
      call ipgeti ('eirene_use_recyceir', use_recyceir)
c*** Initialise tracing routines
      if(ank_tracing.eq.0) return
      if(filedata(ntrc_bln)) call b2blnci(nx,ny,ns,dtim,iline)
      if(filedata(ntrc_bls)) call b2blnsi(nx,ny,ns,dtim,iline)
      if(filedata(ntrc_ble)) call b2blnei(nx,ny,ns,BoRiS,iline)
      if(filedata(ntrc_blm)) call b2blnmi(nx,ny,ns,iline)
      if(filedata(ntrc_int)) call b2wrinti(nx,ny,ns,nnreg(0),iline)
      if(filedata(ntrc_srt)) call b2wrsrti(nx,ny,ns,nnreg(0),iline)
      if(filedata(ntrc_sep)) call b2wrsepi(nx,ny,ns,BoRiS,iline)
      if(filedata(ntrc_usr))
     > call b2usrtri(nx,ny,ns,nxtl,nxtr,BoRiS,iline)
      tracing_initialised = .true.
      return
c======================================================================
      entry b2trcc(nput,nx,ny,ns,nxtl,nxtr,ismain,dtim,BoRiS)
c*** Truncate to match savefile
      if(filedata(ntrc_bln)) then
        call b2blncc(nx,ny,ns,dtim,iline)
        if (iline.gt.0)
     >   write(nput,'(a,i8)') trim(fn_ntrc(ntrc_bln)), iline
      end if
      if(filedata(ntrc_bls)) then
        call b2blnsc(nx,ny,ns,dtim,iline)
        if (iline.gt.0)
     >   write(nput,'(a,i8)') trim(fn_ntrc(ntrc_bls)), iline
      end if
      if(filedata(ntrc_ble)) then
        call b2blnec(nx,ny,ns,BoRiS,iline)
        if (iline.gt.0)
     >  write(nput,'(a,i8)') trim(fn_ntrc(ntrc_ble)), iline
      end if
      if(filedata(ntrc_blm)) then
        call b2blnmc(nx,ny,ns,iline)
        if (iline.gt.0)
     >  write(nput,'(a,i8)') trim(fn_ntrc(ntrc_blm)), iline
      end if
      if(filedata(ntrc_int)) then
        call b2wrintc(nx,ny,ns,nnreg(0),iline)
        if (iline.gt.0)
     >  write(nput,'(a,i8)') trim(fn_ntrc(ntrc_int)), iline
      end if
      if(filedata(ntrc_srt)) then
        call b2wrsrtc(nx,ny,ns,nnreg(0),iline)
        if (iline.gt.0)
     >  write(nput,'(a,i8)') trim(fn_ntrc(ntrc_srt)), iline
      end if
      if(filedata(ntrc_sep)) then
        call b2wrsepc(nx,ny,ns,BoRiS,iline)
        if (iline.gt.0)
     >  write(nput,'(a,i8)') trim(fn_ntrc(ntrc_sep)), iline
      end if
      if(filedata(ntrc_usr)) then
        call b2usrtrcc(nx,ny,ns,nxtl,nxtr,BoRiS,iline)
        if (iline.gt.0)
     >  write(nput,'(a,i8)') trim(fn_ntrc(ntrc_usr)), iline
      end if
      return
c======================================================================
      entry b2trcs(nput,nx,ny,ns,nxtl,nxtr,ismain,dtim,BoRiS)
c*** - to be called BEFORE each timestep
      if(ank_tracing.eq.0) return
      if(filedata(ntrc_bln)) call b2blncs(nx,ny,ns,dtim,iline)
      if(filedata(ntrc_bls)) call b2blnss(nx,ny,ns,dtim,iline)
      if(filedata(ntrc_ble)) call b2blnes(nx,ny,ns,BoRiS,iline)
      if(filedata(ntrc_blm)) call b2blnms(nx,ny,ns,iline)
      if(filedata(ntrc_int)) call b2wrints(nx,ny,ns,nnreg(0),iline)
      if(filedata(ntrc_srt)) call b2wrsrts(nx,ny,ns,nnreg(0),iline)
      if(filedata(ntrc_sep)) call b2wrseps(nx,ny,ns,BoRiS,iline)
      if(filedata(ntrc_usr))
     > call b2usrtrs(nx,ny,ns,nxtl,nxtr,BoRiS,iline)
      return
c======================================================================
      entry b2trcf(nput,nx,ny,ns,nxtl,nxtr,ismain,dtim,BoRiS)
c*** Finalization call
      if (.not.tracing_initialised) return
      call dealloc_b2mod_user
      if(ank_tracing.eq.0) return
      if(filedata(ntrc_bln)) call b2blncf(nx,ny,ns,dtim,iline)
      if(filedata(ntrc_bls)) call b2blnsf(nx,ny,ns,dtim,iline)
      if(filedata(ntrc_ble)) call b2blnef(nx,ny,ns,BoRiS,iline)
      if(filedata(ntrc_blm)) call b2blnmf(nx,ny,ns,iline)
      if(filedata(ntrc_int)) call b2wrintf(nx,ny,ns,nnreg(0),iline)
      if(filedata(ntrc_srt)) call b2wrsrtf(nx,ny,ns,nnreg(0),iline)
      if(filedata(ntrc_sep)) call b2wrsepf(nx,ny,ns,BoRiS,iline)
      if(filedata(ntrc_usr))
     > call b2usrtrf(nx,ny,ns,nxtl,nxtr,BoRiS,iline)
      tracing_initialised = .false.
      return
c======================================================================
      end subroutine b2trace
