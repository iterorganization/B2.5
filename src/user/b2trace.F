      subroutine b2trace(nx,ny,ns,ismain,dtim,BoRiS)
      use b2mod_types
      use b2mod_diag
      use b2mod_indirect
      use b2mod_user_namelist
      use b2mod_neutrals_namelist
csw 07feb2011
      use eirmod_extrab25
csw
c
c  version : 29.08.2000 00:26
c
c======================================================================
c*** Output of tracing information at the end of each outer iteration
c----------------------------------------------------------------------
      implicit none
c
csw 07feb2011
csw#     include <COUPLE/EIRDIAG.F>
c
      integer nx, ny, ns, ismain
      real (kind=R8) :: dtim, BoRiS
      integer ix, iy, is, istra
      logical, save :: tracing_initialised
      data tracing_initialised /.false./
c----------------------------------------------------------------------
c
c*** Correct the neutral fluxes according to recyc values
c
      write(*,"(a,5i4)") ' nspecies,nnatmi,nnstsi,nnlimi,nnmoli=',   !srv 18.05.09
     1 nspecies,nnatmi,nnstsi,nnlimi,nnmoli
      if(nspecies.eq.nnatmi) then
        do ix=1,nnstsi+nnlimi
csw        do ix=1,nlmpgs
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
              wldnek(ix,0)=wldnek(ix,0)+
     +         wldnek(ix,istra)*recyc(ismain,istra)
              wldnep(ix,0)=wldnep(ix,0)+
     +         wldnep(ix,istra)*recyc(ismain,istra)
              iy=ions_list(ismain)
              do is=1,nnatmi
                wldna(ix,is,0)=wldna(ix,is,0)+
     +           wldna(ix,is,istra)*recyc(iy,istra)
                wldra(ix,is,0)=wldra(ix,is,0)+
     +           wldra(ix,is,istra)*recyc(iy,istra)
                ewlda(ix,is,0)=ewlda(ix,is,0)+
     +           ewlda(ix,is,istra)*recyc(iy,istra)
                iy=iy+nfluids(is)
              end do
              do is=1,nnmoli
                iy=eb2spcr(lmolscl(is))+1
                wldnm(ix,is,0)=wldnm(ix,is,0)+
     +           wldnm(ix,is,istra)*recyc(iy,istra)
                wldrm(ix,is,0)=wldrm(ix,is,0)+
     +           wldrm(ix,is,istra)*recyc(iy,istra)
                ewldm(ix,is,0)=ewldm(ix,is,0)+
     +           ewldm(ix,is,istra)*recyc(iy,istra)
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
c*** Momentum balance
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
      end
