      module b2mod_trace
      use b2mod_types
      use b2mod_diag
      use b2mod_blnc
      use b2mod_blne
      use b2mod_blnm
      use b2mod_blns
      use b2mod_wrint
      use b2mod_wrsep
      use b2mod_wrsrt
      use b2mod_usrtrc
      use b2mod_usrspb
      implicit none

      private

      public :: b2trace, b2trci, b2trcf, b2trcs, b2trcc

      integer, save :: use_recyceir
      data use_recyceir /1/
      logical, save :: tracing_initialised
      data tracing_initialised /.false./

      contains

      subroutine b2trace(nx,ny,ns,nxtl,nxtr,nregionv,ismain,dtim,BoRiS)
      use b2mod_b2cmpa
      use b2mod_indirect
      use b2mod_user_namelist
      use b2mod_neutrals_namelist
#ifdef B25_EIRENE
      use eirmod_wneutrals
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
      integer nx, ny, ns, nxtl, nxtr, nregionv, ismain
      real (kind=R8) :: dtim, BoRiS
#ifdef B25_EIRENE
      integer ix, iy, is, istra
      real (kind=R8) :: recyctmp
c----------------------------------------------------------------------
c
c*** Correct the neutral fluxes according to recyc values
c
      write(*,'(a,6i5)')
     , 'nspecies,nnatmi,nnstsi,nnlimi,nnmoli,nnplsi =',
     1  nspecies,nnatmi,nnstsi,nnlimi,nnmoli,nnplsi
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
      if(filedata(ntrc_bln)) call b2blnc(nx,ny,ns,dtim)
c*** Particle balance in SPb style (for use in the Method of effective source for code speed-up,
c***  see E. Kaveeva et al 2018 Nucl. Fusion 58 126018)
      if(filedata(ntrc_bls)) call b2blns(nx,ny,ns,dtim)
c*** Power balance
      if(filedata(ntrc_ble)) call b2blne(nx,ny,ns,BoRiS)
c*** Parallel momentum balance
      if(filedata(ntrc_blm)) call b2blnm(nx,ny,ns)
c*** Additional integrals
      if(filedata(ntrc_int)) call b2wrint(nx,ny,ns,nregionv)
      if(filedata(ntrc_srt)) call b2wrsrt(nx,ny,ns)
c*** Separatrix data
      if(filedata(ntrc_sep)) call b2wrsep(nx,ny,ns,BoRiS)
c*** User-specific data
      if(filedata(ntrc_usr)) call b2usrtrc(nx,ny,ns,nxtl,nxtr,BoRiS)
c*** User-specific data
      if(filedata(ntrc_uss)) call b2usstrc(nx,ny,ns,BoRiS)       ! IYS 22.03.2022
c
      return
      end subroutine b2trace
c======================================================================
      subroutine b2trci(nx,ny,ns,nxtl,nxtr,nregionv)
      implicit none
      integer, intent(in) :: nx, ny, ns, nxtl,nxtr, nregionv
c*** Initialisation call
      if (tracing_initialised) return
      call alloc_b2mod_user(nx,ns,nlim,nmol,ntns)
      call ipgeti ('eirene_use_recyceir', use_recyceir)
c*** Initialise tracing routines
      if(ank_tracing.eq.0) return
      if(filedata(ntrc_bln)) call b2blnci
      if(filedata(ntrc_bls)) call b2blnsi
      if(filedata(ntrc_ble)) call b2blnei(ny)
      if(filedata(ntrc_blm)) call b2blnmi
      if(filedata(ntrc_int)) call b2wrinti(nregionv)
      if(filedata(ntrc_srt)) call b2wrsrti(nregionv)
      if(filedata(ntrc_sep)) call b2wrsepi(nx,ny,ns)
      if(filedata(ntrc_usr)) call b2usrtri(nx,ny,ns)
      if(filedata(ntrc_uss)) call b2usstri(nx,ny,ns,nxtl,nxtr) ! IYS 22.03.2022
      tracing_initialised = .true.
      return
      end subroutine b2trci
c======================================================================
      subroutine b2trcc(nput,nx,nxtl,nxtr)
      implicit none
      integer, intent(in) :: nput, nx, nxtl, nxtr
      integer iline
c*** Truncate to match savefile
      if(filedata(ntrc_bln)) then
        call b2blncc(iline)
        if (iline.gt.0)
     >   write(nput,'(a,i8)') trim(fn_ntrc(ntrc_bln)), iline
      end if
      if(filedata(ntrc_bls)) then
        call b2blnsc(iline)
        if (iline.gt.0)
     >   write(nput,'(a,i8)') trim(fn_ntrc(ntrc_bls)), iline
      end if
      if(filedata(ntrc_ble)) then
        call b2blnec(iline)
        if (iline.gt.0)
     >  write(nput,'(a,i8)') trim(fn_ntrc(ntrc_ble)), iline
      end if
      if(filedata(ntrc_blm)) then
        call b2blnmc(iline)
        if (iline.gt.0)
     >  write(nput,'(a,i8)') trim(fn_ntrc(ntrc_blm)), iline
      end if
      if(filedata(ntrc_int)) then
        call b2wrintc(iline)
        if (iline.gt.0)
     >  write(nput,'(a,i8)') trim(fn_ntrc(ntrc_int)), iline
      end if
      if(filedata(ntrc_srt)) then
        call b2wrsrtc(iline)
        if (iline.gt.0)
     >  write(nput,'(a,i8)') trim(fn_ntrc(ntrc_srt)), iline
      end if
      if(filedata(ntrc_sep)) then
        call b2wrsepc(iline)
        if (iline.gt.0)
     >  write(nput,'(a,i8)') trim(fn_ntrc(ntrc_sep)), iline
      end if
      if(filedata(ntrc_usr)) then
        call b2usrtrcc(nx,nxtl,nxtr,iline)
        if (iline.gt.0)
     >  write(nput,'(a,i8)') trim(fn_ntrc(ntrc_usr)), iline
      end if
      if(filedata(ntrc_uss)) then      ! IYS 22.03.2022
        call b2usstrcc(iline)
        if (iline.gt.0)
     >  write(nput,'(a,i8)') trim(fn_ntrc(ntrc_uss)), iline
      end if
      return
      end subroutine b2trcc
c======================================================================
      subroutine b2trcs
      implicit none
c*** - to be called BEFORE each timestep
      if(ank_tracing.eq.0) return
      if(filedata(ntrc_bln)) call b2blncs
      if(filedata(ntrc_bls)) call b2blnss
      if(filedata(ntrc_ble)) call b2blnes
      if(filedata(ntrc_blm)) call b2blnms
      if(filedata(ntrc_int)) call b2wrints
      if(filedata(ntrc_srt)) call b2wrsrts
      if(filedata(ntrc_sep)) call b2wrseps
      if(filedata(ntrc_usr)) call b2usrtrs
      if(filedata(ntrc_uss)) call b2usstrs      ! IYS 22.03.2022
      return
      end subroutine b2trcs
c======================================================================
      subroutine b2trcf
      implicit none
c*** Finalization call
      if (.not.tracing_initialised) return
      call dealloc_b2mod_user
      if(ank_tracing.eq.0) return
      if(filedata(ntrc_bln)) call b2blncf
      if(filedata(ntrc_bls)) call b2blnsf
      if(filedata(ntrc_ble)) call b2blnef
      if(filedata(ntrc_blm)) call b2blnmf
      if(filedata(ntrc_int)) call b2wrintf
      if(filedata(ntrc_srt)) call b2wrsrtf
      if(filedata(ntrc_sep)) call b2wrsepf
      if(filedata(ntrc_usr)) call b2usrtrf
      if(filedata(ntrc_uss)) call b2usstrf      ! IYS 22.03.2022
      tracing_initialised = .false.
      return
      end subroutine b2trcf
c======================================================================
      end module b2mod_trace
