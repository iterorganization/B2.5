      module b2mod_sources
      use b2mod_types
      use b2mod_version
      implicit none
      real (kind=R8), allocatable, save ::
     *  sch(:,:,:),
     *  she(:,:,:),
     *  shi(:,:,:),
     *  sne(:,:,:),
     *  smo(:,:,:,:),
     *  shef(:,:,:),                                           !srv 10.11.02                
     *  shif(:,:,:),                                           !srv 10.11.02 
     *  smof(:,:,:,:),                                         !srv 10.11.02 
     *  smq(:,:,:,:),
     *  sna(:,:,:,:),
     *  rsana(:,:,:),
     *  rsahi(:,:,:),
     *  rsamo(:,:,:),
     *  rrana(:,:,:),
     *  rrahi(:,:,:),
     *  rramo(:,:,:),
     *  rqahe(:,:,:),
     *  rqrad(:,:,:),
     *  rqbrm(:,:,:),
     *  rcxna(:,:,:),
     *  rcxhi(:,:,:),
     *  rcxmo(:,:,:)
      real (kind=R8), allocatable, save ::
     *  b2stbr_sna(:,:,:),
     *  b2stbr_smo(:,:,:),
     *  b2stbr_sne(:,:),
     *  b2stbr_she(:,:),
     *  b2stbr_shi(:,:),
     *  b2stbr_sch(:,:),
     *  b2stbc_sna(:,:,:),
     *  b2stbc_smo(:,:,:),
     *  b2stbc_sne(:,:),
     *  b2stbc_she(:,:),
     *  b2stbc_shi(:,:),
     *  b2stbc_sch(:,:),
     *  b2sihs_divue(:,:),
     *  b2sihs_divua(:,:),
     *  b2sihs_exbe(:,:),
     *  b2sihs_exba(:,:),
     *  b2sihs_visa(:,:),
     *  b2sihs_joule(:,:),
     *  b2sihs_fraa(:,:),
     *  b2sihs_str(:,:),
     1  b2npmo_smaf(:,:,:,:),
     1  b2npmo_smag(:,:,:,:),
     1  b2npmo_smav(:,:,:,:),
     *  ext_sna(:,:,:),
     *  ext_smo(:,:,:),
     *  ext_sne(:,:),
     *  ext_she(:,:),
     *  ext_shi(:,:),
     *  ext_sch(:,:)
      real (kind=R8), allocatable, save ::                             !srv 11.09.09 16.07.10
     *  snasm(:,:,:,:), smosm(:,:,:,:), 
     &  shesm(:,:,:), shism(:,:,:),
     *  snadt(:,:,:,:), smodt(:,:,:,:), 
     &  shedt(:,:,:), shidt(:,:,:),
     &  schdt(:,:,:)

      contains

      subroutine alloc_b2mod_sources(nxd,nyd,nsd)
      implicit none
      integer nxd,nyd,nsd

      allocate(sch(-1:nxd,-1:nyd,0:3))
      allocate(she(-1:nxd,-1:nyd,0:3))
      allocate(sne(-1:nxd,-1:nyd,0:1))
      allocate(shi(-1:nxd,-1:nyd,0:3))
      allocate(smo(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(sna(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(shef(-1:nxd,-1:nyd,0:3))                        !sv 10.11.02
      allocate(shif(-1:nxd,-1:nyd,0:3))                        !sv 10.11.02
      allocate(smof(-1:nxd,-1:nyd,0:3,0:nsd-1))                !sv 10.11.02
      allocate(smq(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2npmo_smag(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2npmo_smaf(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2npmo_smav(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(rsana(-1:nxd,-1:nyd,0:nsd-1))
      allocate(rsahi(-1:nxd,-1:nyd,0:nsd-1))
      allocate(rsamo(-1:nxd,-1:nyd,0:nsd-1))
      allocate(rrana(-1:nxd,-1:nyd,0:nsd-1))
      allocate(rrahi(-1:nxd,-1:nyd,0:nsd-1))
      allocate(rramo(-1:nxd,-1:nyd,0:nsd-1))
      allocate(rqahe(-1:nxd,-1:nyd,0:nsd-1))
      allocate(rqrad(-1:nxd,-1:nyd,0:nsd-1))
      allocate(rqbrm(-1:nxd,-1:nyd,0:nsd-1))
      allocate(rcxna(-1:nxd,-1:nyd,0:nsd-1))
      allocate(rcxhi(-1:nxd,-1:nyd,0:nsd-1))
      allocate(rcxmo(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbr_sna(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbr_smo(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbc_sne(-1:nxd,-1:nyd))
      allocate(b2stbr_she(-1:nxd,-1:nyd))
      allocate(b2stbr_shi(-1:nxd,-1:nyd))
      allocate(b2stbr_sch(-1:nxd,-1:nyd))
      allocate(b2stbc_sna(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbc_smo(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbr_sne(-1:nxd,-1:nyd))
      allocate(b2stbc_she(-1:nxd,-1:nyd))
      allocate(b2stbc_shi(-1:nxd,-1:nyd))
      allocate(b2stbc_sch(-1:nxd,-1:nyd))
      allocate(b2sihs_divue(-1:nxd,-1:nyd))
      allocate(b2sihs_divua(-1:nxd,-1:nyd))
      allocate(b2sihs_exbe(-1:nxd,-1:nyd))
      allocate(b2sihs_exba(-1:nxd,-1:nyd))
      allocate(b2sihs_visa(-1:nxd,-1:nyd))
      allocate(b2sihs_joule(-1:nxd,-1:nyd))
      allocate(b2sihs_fraa(-1:nxd,-1:nyd))
      allocate(b2sihs_str(-1:nxd,-1:nyd))
      allocate(ext_sna(-1:nxd,-1:nyd,0:nsd-1))
      allocate(ext_smo(-1:nxd,-1:nyd,0:nsd-1))
      allocate(ext_sne(-1:nxd,-1:nyd))
      allocate(ext_she(-1:nxd,-1:nyd))
      allocate(ext_shi(-1:nxd,-1:nyd))
      allocate(ext_sch(-1:nxd,-1:nyd))
      allocate(schdt(-1:nxd,-1:nyd,0:3))                        !srv 16.07.10 {
      allocate(shedt(-1:nxd,-1:nyd,0:3))
      allocate(shidt(-1:nxd,-1:nyd,0:3))
      allocate(smodt(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(snadt(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(shesm(-1:nxd,-1:nyd,0:3))
      allocate(shism(-1:nxd,-1:nyd,0:3))
      allocate(smosm(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(snasm(-1:nxd,-1:nyd,0:1,0:nsd-1))                 !srv 16.07.10 }

      rsana=0.0_R8
      rsahi=0.0_R8
      rsamo=0.0_R8
      rrana=0.0_R8
      rrahi=0.0_R8
      rramo=0.0_R8
      rqahe=0.0_R8
      rqrad=0.0_R8
      rqbrm=0.0_R8
      rcxna=0.0_R8
      rcxhi=0.0_R8
      rcxmo=0.0_R8
      ext_sna=0.0_R8
      ext_smo=0.0_R8
      ext_sne=0.0_R8
      ext_she=0.0_R8
      ext_shi=0.0_R8
      ext_sch=0.0_R8

      return
      end subroutine alloc_b2mod_sources

      subroutine dealloc_b2mod_sources
      implicit none

      deallocate(sch)
      deallocate(she)
      deallocate(sne)
      deallocate(shi)
      deallocate(smo)
      deallocate(shef)                        !sv 10.11.02
      deallocate(shif)                        !sv 10.11.02
      deallocate(smof)                        !sv 10.11.02
      deallocate(smq)
      deallocate(b2npmo_smag)
      deallocate(b2npmo_smaf)
      deallocate(b2npmo_smav)
      deallocate(sna)
      deallocate(rsana)
      deallocate(rsahi)
      deallocate(rsamo)
      deallocate(rrana)
      deallocate(rrahi)
      deallocate(rramo)
      deallocate(rqahe)
      deallocate(rqrad)
      deallocate(rqbrm)
      deallocate(rcxna)
      deallocate(rcxhi)
      deallocate(rcxmo)
      deallocate(b2stbr_sna)
      deallocate(b2stbr_smo)
      deallocate(b2stbr_sne)
      deallocate(b2stbr_she)
      deallocate(b2stbr_shi)
      deallocate(b2stbr_sch)
      deallocate(b2stbc_sna)
      deallocate(b2stbc_smo)
      deallocate(b2stbc_sne)
      deallocate(b2stbc_she)
      deallocate(b2stbc_shi)
      deallocate(b2stbc_sch)
      deallocate(b2sihs_divue)
      deallocate(b2sihs_divua)
      deallocate(b2sihs_exbe)
      deallocate(b2sihs_exba)
      deallocate(b2sihs_visa)
      deallocate(b2sihs_joule)
      deallocate(b2sihs_fraa)
      deallocate(b2sihs_str)
      deallocate(ext_sna)
      deallocate(ext_sne)
      deallocate(ext_smo)
      deallocate(ext_she)
      deallocate(ext_shi)
      deallocate(ext_sch)
      deallocate(schdt)                        !srv 16.07.10 {
      deallocate(shedt)
      deallocate(shidt)
      deallocate(smodt)
      deallocate(snadt)
      deallocate(shesm)
      deallocate(shism)
      deallocate(smosm)
      deallocate(snasm)                         !srv 16.07.10 }

      return
      end subroutine dealloc_b2mod_sources

      subroutine read_b2mod_sources(nread)
      implicit none
      integer nread

      call cfrure (nread, size(sch), sch, 'sch')
      call cfrure (nread, size(she), she, 'she')
      call cfrure (nread, size(shi), shi, 'shi')
      call cfrure (nread, size(smo), smo, 'smo')
      call cfrure (nread, size(smq), smq, 'smq')
      call cfrure (nread, size(sna), sna, 'sna')
      if(b2fplasma_version.gt.'00.000.000') then
        call cfrure (nread, size(rsana), rsana, 'rsana')
        call cfrure (nread, size(rsahi), rsahi, 'rsahi')
        if(b2fplasma_version.gt.'01.001.003') then
          call cfrure (nread, size(rsamo), rsamo, 'rsamo')
        endif
        call cfrure (nread, size(rrana), rrana, 'rrana')
        call cfrure (nread, size(rrahi), rrahi, 'rrahi')
        if(b2fplasma_version.gt.'01.001.003') then
          call cfrure (nread, size(rramo), rramo, 'rramo')
        endif
        call cfrure (nread, size(rqahe), rqahe, 'rqahe')
        if(b2fplasma_version.ge.'01.001.016') then
          call cfrure (nread, size(rqrad), rqrad, 'rqrad')
        endif
        if(b2fplasma_version.ge.'01.001.019') then
          call cfrure (nread, size(rqbrm), rqbrm, 'rqbrm')
        endif
        call cfrure (nread, size(rcxna), rcxna, 'rcxna')
        call cfrure (nread, size(rcxhi), rcxhi, 'rcxhi')
        if(b2fplasma_version.gt.'01.001.003') then
          call cfrure (nread, size(rcxmo), rcxmo, 'rcxmo')
        endif
      endif
      if(b2fplasma_version.gt.'01.001.000') then
        call cfrure (nread, size(b2stbr_sna), b2stbr_sna, 'b2stbr_sna')
        call cfrure (nread, size(b2stbr_smo), b2stbr_smo, 'b2stbr_smo')
        call cfrure (nread, size(b2stbr_she), b2stbr_she, 'b2stbr_she')
        call cfrure (nread, size(b2stbr_shi), b2stbr_shi, 'b2stbr_shi')
        call cfrure (nread, size(b2stbr_sch), b2stbr_sch, 'b2stbr_sch')
      endif
      if(b2fplasma_version.gt.'01.001.010') then
        call cfrure (nread, size(b2stbc_sna), b2stbc_sna, 'b2stbc_sna')
        call cfrure (nread, size(b2stbc_smo), b2stbc_smo, 'b2stbc_smo')
        call cfrure (nread, size(b2stbc_she), b2stbc_she, 'b2stbc_she')
        call cfrure (nread, size(b2stbc_shi), b2stbc_shi, 'b2stbc_shi')
        call cfrure (nread, size(b2stbc_sch), b2stbc_sch, 'b2stbc_sch')
      endif
      if(b2fplasma_version.gt.'01.001.002') then
        call cfrure (nread, size(b2sihs_divue), b2sihs_divue,
     1   'b2sihs_divue') 
        call cfrure (nread, size(b2sihs_divua), b2sihs_divua,
     1   'b2sihs_divua') 
        call cfrure (nread, size(b2sihs_exbe),  b2sihs_exbe,
     1   'b2sihs_exbe')  
        call cfrure (nread, size(b2sihs_exba),  b2sihs_exba,
     1   'b2sihs_exba')  
        call cfrure (nread, size(b2sihs_visa),  b2sihs_visa,
     1   'b2sihs_visa')  
        call cfrure (nread, size(b2sihs_joule), b2sihs_joule,
     1   'b2sihs_joule') 
        call cfrure (nread, size(b2sihs_fraa),  b2sihs_fraa,
     1   'b2sihs_fraa')  
        call cfrure (nread, size(b2sihs_str),   b2sihs_str,
     1   'b2sihs_str')   
      endif
      if(b2fplasma_version.gt.'01.001.005') then
        call cfrure (nread, size(b2npmo_smaf),  b2npmo_smaf,
     1   'b2npmo_smaf')
        call cfrure (nread, size(b2npmo_smag),  b2npmo_smag,
     1   'b2npmo_smag')
        call cfrure (nread, size(b2npmo_smav),  b2npmo_smav,
     1   'b2npmo_smav')
      endif
      if(b2fplasma_version.gt.'01.001.012') then
        call cfrure (nread, size(ext_sna), ext_sna, 'ext_sna')
        call cfrure (nread, size(ext_smo), ext_smo, 'ext_smo')
        call cfrure (nread, size(ext_she), ext_she, 'ext_she')
        call cfrure (nread, size(ext_shi), ext_shi, 'ext_shi')
        call cfrure (nread, size(ext_sch), ext_sch, 'ext_sch')
      endif

      return
      end subroutine read_b2mod_sources

      subroutine write_b2mod_sources(nwrite)
      implicit none
      integer nwrite

      call cfwure (nwrite, size(sch), sch, 'sch')
      call cfwure (nwrite, size(she), she, 'she')
      call cfwure (nwrite, size(shi), shi, 'shi')
      call cfwure (nwrite, size(smo), smo, 'smo')
      call cfwure (nwrite, size(smq), smq, 'smq')
      call cfwure (nwrite, size(sna), sna, 'sna')
      if(newversion.gt.'00.000.000') then
        write(*,*) newversion,'rsana --- rcxhi'
        call cfwure (nwrite, size(rsana), rsana, 'rsana')
        call cfwure (nwrite, size(rsahi), rsahi, 'rsahi')
        if(newversion.gt.'01.001.003') then
          write(*,*) newversion,'rsamo --- rcxmo'
          call cfwure (nwrite, size(rsamo), rsamo, 'rsamo')
        endif
        call cfwure (nwrite, size(rrana), rrana, 'rrana')
        call cfwure (nwrite, size(rrahi), rrahi, 'rrahi')
        if(newversion.gt.'01.001.003') then
          call cfwure (nwrite, size(rramo), rramo, 'rramo')
        endif
        call cfwure (nwrite, size(rqahe), rqahe, 'rqahe')
        if(newversion.ge.'01.001.016') then
          call cfwure (nwrite, size(rqrad), rqrad, 'rqrad')
        endif
        if(newversion.ge.'01.001.019') then
          call cfwure (nwrite, size(rqbrm), rqbrm, 'rqbrm')
        endif
        call cfwure (nwrite, size(rcxna), rcxna, 'rcxna')
        call cfwure (nwrite, size(rcxhi), rcxhi, 'rcxhi')
        if(newversion.gt.'01.001.003') then
          call cfwure (nwrite, size(rcxmo), rcxmo, 'rcxmo')
        endif
      endif
      if(newversion.gt.'01.001.000') then
        write(*,*) newversion,'b2stbr_sna --- _sch'
        call cfwure (nwrite, size(b2stbr_sna), b2stbr_sna, 'b2stbr_sna')
        call cfwure (nwrite, size(b2stbr_smo), b2stbr_smo, 'b2stbr_smo')
        call cfwure (nwrite, size(b2stbr_she), b2stbr_she, 'b2stbr_she')
        call cfwure (nwrite, size(b2stbr_shi), b2stbr_shi, 'b2stbr_shi')
        call cfwure (nwrite, size(b2stbr_sch), b2stbr_sch, 'b2stbr_sch')
      endif
      if(newversion.gt.'01.001.010') then
        write(*,*) newversion,'b2stbc_sna --- _sch'
        call cfwure (nwrite, size(b2stbc_sna), b2stbc_sna, 'b2stbc_sna')
        call cfwure (nwrite, size(b2stbc_smo), b2stbc_smo, 'b2stbc_smo')
        call cfwure (nwrite, size(b2stbc_she), b2stbc_she, 'b2stbc_she')
        call cfwure (nwrite, size(b2stbc_shi), b2stbc_shi, 'b2stbc_shi')
        call cfwure (nwrite, size(b2stbc_sch), b2stbc_sch, 'b2stbc_sch')
      endif
      if(newversion.gt.'01.001.002') then
        write(*,*) newversion,'b2sihs_divue --- _str'
        call cfwure (nwrite, size(b2sihs_divue),b2sihs_divue,
     1   'b2sihs_divue') 
        call cfwure (nwrite, size(b2sihs_divua),b2sihs_divua,
     1   'b2sihs_divua') 
        call cfwure (nwrite, size(b2sihs_exbe), b2sihs_exbe,
     1   'b2sihs_exbe')  
        call cfwure (nwrite, size(b2sihs_exba), b2sihs_exba,
     1   'b2sihs_exba')  
        call cfwure (nwrite, size(b2sihs_visa), b2sihs_visa,
     1   'b2sihs_visa')  
        call cfwure (nwrite, size(b2sihs_joule),b2sihs_joule,
     1   'b2sihs_joule') 
        call cfwure (nwrite, size(b2sihs_fraa), b2sihs_fraa,
     1   'b2sihs_fraa')  
        call cfwure (nwrite, size(b2sihs_str),  b2sihs_str,
     1   'b2sihs_str')   
      endif
      if(newversion.gt.'01.001.005') then
        call cfwure (nwrite, size(b2npmo_smaf),  b2npmo_smaf,
     1   'b2npmo_smaf')
        call cfwure (nwrite, size(b2npmo_smag),  b2npmo_smag,
     1   'b2npmo_smag')
        call cfwure (nwrite, size(b2npmo_smav),  b2npmo_smav,
     1   'b2npmo_smav')
      endif
      if(newversion.gt.'01.001.012') then
        call cfwure (nwrite, size(ext_sna), ext_sna, 'ext_sna')
        call cfwure (nwrite, size(ext_smo), ext_smo, 'ext_smo')
        call cfwure (nwrite, size(ext_she), ext_she, 'ext_she')
        call cfwure (nwrite, size(ext_shi), ext_shi, 'ext_shi')
        call cfwure (nwrite, size(ext_sch), ext_sch, 'ext_sch')
      endif

      return
      end subroutine write_b2mod_sources

      end module b2mod_sources
