      subroutine b2trace(nx,ny,ns,ismain,dtim,BoRiS)
      use b2mod_types
      use b2mod_diag
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
c----------------------------------------------------------------------
      implicit none
      integer nx, ny, ns, ismain
      real (kind=R8) :: dtim, BoRiS
      integer ix, iy, is, istra
      integer, save :: use_recyceir
      real (kind=R8) :: recyctmp
      logical, save :: tracing_initialised
      data tracing_initialised /.false./
      data use_recyceir /0/
c----------------------------------------------------------------------
c
c*** Correct the neutral fluxes according to recyc values
c
      write(*,'(a,6i5)') 'nspecies,nnatmi,nnstsi,nnlimi,nnmoli,nnplsi=',
     1 nspecies,nnatmi,nnstsi,nnlimi,nnmoli,nnplsi
      if(nspecies.eq.nnatmi) then
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
     +       (ix.gt.nnlimi.and.ix.le.nnlimi+nnstsi)) then
            do istra=1,nstrai
              if(use_recyceir.gt.0) then
                recyctmp=recyceir(istra)
              else
                recyctmp=recyc(ismain,istra)
              endif
              wldnek(ix,0)=wldnek(ix,0)+
     +         wldnek(ix,istra)*recyctmp
              wldnep(ix,0)=wldnep(ix,0)+
     +         wldnep(ix,istra)*recyctmp
              do is=1,nnatmi
                iy=b2espcr(is)
                if(use_recyceir.eq.0) recyctmp=recyc(iy,istra)
                wldna(ix,is,0)=wldna(ix,is,0)+
     +           wldna(ix,is,istra)*recyctmp
                wldra(ix,is,0)=wldra(ix,is,0)+
     +           wldra(ix,is,istra)*recyctmp
                ewlda(ix,is,0)=ewlda(ix,is,0)+
     +           ewlda(ix,is,istra)*recyctmp
                iy=iy+nfluids(is)
              end do
              do is=1,nnmoli
                if (lmolscl(is).gt.0 .and. lmolscl(is).le.nnatmi) then
                  iy=eb2spcr(lmolscl(is))+1
                  if(use_recyceir.eq.0) recyctmp=recyc(iy,istra)
                  wldnm(ix,is,0)=wldnm(ix,is,0)+
     +             wldnm(ix,is,istra)*recyctmp
                  wldrm(ix,is,0)=wldrm(ix,is,0)+
     +             wldrm(ix,is,istra)*recyctmp
                  ewldm(ix,is,0)=ewldm(ix,is,0)+
     +             ewldm(ix,is,istra)*recyctmp
                else if (use_recyceir.gt.0) then !here we are rescaling by stratum, not by species
                  wldnm(ix,is,0)=wldnm(ix,is,0)+
     +             wldnm(ix,is,istra)*recyctmp
                  wldrm(ix,is,0)=wldrm(ix,is,0)+
     +             wldrm(ix,is,istra)*recyctmp
                  ewldm(ix,is,0)=ewldm(ix,is,0)+
     +             ewldm(ix,is,istra)*recyctmp
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
      else
        write (*,"(a,5i4)")
     &  'b2trace: no source rescaling: natmi,nspecies=', nnatmi,nspecies  !srv 18.05.09
      end if
c*** Particle balance
      call b2blnc(nx,ny,ns,dtim)
c*** Power balance
      call b2blne(nx,ny,ns,BoRiS)
c*** Parallel momentum balance
      call b2blnm(nx,ny,ns)
c*** Additional integrals
      call b2wrint(nx,ny,ns,nnreg(0))
c*** Separatrix data
      call b2wrsep(nx,ny,ns,BoRiS)
c*** User-specific data
      call b2usrtrc(nx,ny,ns,BoRiS)
c
      return
c======================================================================
      entry b2trci(nx,ny,ns,ismain,dtim,BoRiS)
c*** Initialisation call
      if (tracing_initialised) return
      call alloc_b2mod_user(nlim,nmol)
      call ipgeti ('eirene_use_recyceir', use_recyceir)
      call read_b2mod_user_namelist(nx)
      if(filedata(ntrc_bln)) call b2blnci(nx,ny,ns,dtim)
      if(filedata(ntrc_ble)) call b2blnei(nx,ny,ns,BoRiS)
      if(filedata(ntrc_blm)) call b2blnmi(nx,ny,ns)
      if(filedata(ntrc_int)) call b2wrinti(nx,ny,ns,nnreg(0))
      if(filedata(ntrc_sep)) call b2wrsepi(nx,ny,ns,BoRiS)
      if(filedata(ntrc_usr)) call b2usrtri(nx,ny,ns,BoRiS)
      tracing_initialised = .true.
      return
c======================================================================
      entry b2trcs(nx,ny,ns,ismain,dtim,BoRiS)
c*** - to be called BEFORE each timestep
      call b2blncs(nx,ny,ns,dtim)
      call b2blnes(nx,ny,ns,BoRiS)
      call b2blnms(nx,ny,ns)
      call b2wrints(nx,ny,ns,nnreg(0))
      call b2wrseps(nx,ny,ns,BoRiS)
      call b2usrtrs(nx,ny,ns,BoRiS)
      return
c======================================================================
      entry b2trcf(nx,ny,ns,ismain,dtim,BoRiS)
c*** Finalization call
      if (.not.tracing_initialised) return
      call dealloc_b2mod_user
      if(filedata(ntrc_bln)) call b2blncf(nx,ny,ns,dtim)
      if(filedata(ntrc_ble)) call b2blnef(nx,ny,ns,BoRiS)
      if(filedata(ntrc_blm)) call b2blnmf(nx,ny,ns)
      if(filedata(ntrc_int)) call b2wrintf(nx,ny,ns,nnreg(0))
      if(filedata(ntrc_sep)) call b2wrsepf(nx,ny,ns,BoRiS)
      if(filedata(ntrc_usr)) call b2usrtrf(nx,ny,ns,BoRiS)
      tracing_initialised = .false.
      return
c======================================================================
      end subroutine b2trace
