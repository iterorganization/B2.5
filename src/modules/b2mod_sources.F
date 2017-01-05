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
     *  b2stel_sna_ion(:,:,:), !djm 14/10/16
     *  b2stel_sna_rec(:,:,:),
     *  b2stel_smq_ion(:,:,:),
     *  b2stel_smq_rec(:,:,:),
     *  b2stel_she(:,:),
     *  b2stel_shi_ion(:,:),
     *  b2stel_shi_rec(:,:),
     *  b2stcx_sna(:,:,:),
     *  b2stcx_smq(:,:,:),
     *  b2stcx_shi(:,:),
     *  b2srsm_sna(:,:,:),
     *  b2srsm_smo(:,:,:),
     *  b2srsm_she(:,:),
     *  b2srsm_shi(:,:),
     *  b2srdt_sna(:,:,:),
     *  b2srdt_smo(:,:,:),
     *  b2srdt_she(:,:),
     *  b2srdt_shi(:,:),
     *  b2srdt_sch(:,:),
     *  b2srst_sna(:,:,:),
     *  b2srst_smo(:,:,:),
     *  b2srst_she(:,:),
     *  b2srst_shi(:,:),
     *  b2srst_sch(:,:),
     *  b2stbm_sna(:,:,:),
     *  b2stbm_smo(:,:,:),
     *  b2stbm_sne(:,:),
     *  b2stbm_she(:,:),
     *  b2stbm_shi(:,:),
     *  b2stbm_sch(:,:),
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
     1  smpr(:,:,:),
     1  smpt(:,:,:),
     1  smfr(:,:,:),
     1  smcf(:,:,:),
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
     &  schdt(:,:,:), snedt(:,:,:)

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
      allocate(smpr(-1:nxd,-1:nyd,0:nsd-1))
      allocate(smpt(-1:nxd,-1:nyd,0:nsd-1))
      allocate(smfr(-1:nxd,-1:nyd,0:nsd-1))
      allocate(smcf(-1:nxd,-1:nyd,0:nsd-1))
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
      allocate(b2stbr_sne(-1:nxd,-1:nyd))
      allocate(b2stbr_she(-1:nxd,-1:nyd))
      allocate(b2stbr_shi(-1:nxd,-1:nyd))
      allocate(b2stbr_sch(-1:nxd,-1:nyd))
      allocate(b2stbc_sna(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbc_smo(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbc_sne(-1:nxd,-1:nyd))
      allocate(b2stbc_she(-1:nxd,-1:nyd))
      allocate(b2stbc_shi(-1:nxd,-1:nyd))
      allocate(b2stbc_sch(-1:nxd,-1:nyd))
      allocate(b2stel_sna_ion(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stel_sna_rec(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stel_smq_ion(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stel_smq_rec(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stel_she(-1:nxd,-1:nyd))
      allocate(b2stel_shi_ion(-1:nxd,-1:nyd))
      allocate(b2stel_shi_rec(-1:nxd,-1:nyd))
      allocate(b2stcx_sna(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stcx_smq(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stcx_shi(-1:nxd,-1:nyd))
      allocate(b2srsm_sna(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2srsm_smo(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2srsm_she(-1:nxd,-1:nyd))
      allocate(b2srsm_shi(-1:nxd,-1:nyd))
      allocate(b2srdt_sna(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2srdt_smo(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2srdt_she(-1:nxd,-1:nyd))
      allocate(b2srdt_shi(-1:nxd,-1:nyd))
      allocate(b2srdt_sch(-1:nxd,-1:nyd))     
      allocate(b2srst_sna(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2srst_smo(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2srst_she(-1:nxd,-1:nyd))
      allocate(b2srst_shi(-1:nxd,-1:nyd))
      allocate(b2srst_sch(-1:nxd,-1:nyd))     
      allocate(b2stbm_sna(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbm_smo(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbm_sne(-1:nxd,-1:nyd))
      allocate(b2stbm_she(-1:nxd,-1:nyd))
      allocate(b2stbm_shi(-1:nxd,-1:nyd))
      allocate(b2stbm_sch(-1:nxd,-1:nyd))
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
      allocate(snedt(-1:nxd,-1:nyd,0:1))
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
      b2stbr_sna=0.0_R8
      b2stbr_smo=0.0_R8
      b2stbr_sne=0.0_R8
      b2stbr_she=0.0_R8
      b2stbr_shi=0.0_R8
      b2stbr_sch=0.0_R8
      b2stbc_sna=0.0_R8
      b2stbc_smo=0.0_R8
      b2stbc_sne=0.0_R8
      b2stbc_she=0.0_R8
      b2stbc_shi=0.0_R8
      b2stbc_sch=0.0_R8
      b2stel_sna_ion=0.0_R8
      b2stel_sna_rec=0.0_R8
      b2stel_smq_ion=0.0_R8
      b2stel_smq_rec=0.0_R8
      b2stel_she=0.0_R8
      b2stel_shi_ion=0.0_R8
      b2stel_shi_rec=0.0_R8
      b2stcx_sna=0.0_R8
      b2stcx_smq=0.0_R8
      b2stcx_shi=0.0_R8
      b2srsm_sna=0.0_R8
      b2srsm_smo=0.0_R8
      b2srsm_she=0.0_R8
      b2srsm_shi=0.0_R8
      b2srdt_sna=0.0_R8
      b2srdt_smo=0.0_R8
      b2srdt_she=0.0_R8
      b2srdt_shi=0.0_R8
      b2srdt_sch=0.0_R8
      b2srst_sna=0.0_R8
      b2srst_smo=0.0_R8
      b2srst_she=0.0_R8
      b2srst_shi=0.0_R8
      b2srst_sch=0.0_R8
      b2stbm_sna=0.0_R8
      b2stbm_smo=0.0_R8
      b2stbm_sne=0.0_R8
      b2stbm_she=0.0_R8
      b2stbm_shi=0.0_R8
      b2stbm_sch=0.0_R8
      b2sihs_divue=0.0_R8
      b2sihs_divua=0.0_R8
      b2sihs_exbe=0.0_R8
      b2sihs_exba=0.0_R8
      b2sihs_visa=0.0_R8
      b2sihs_joule=0.0_R8
      b2sihs_fraa=0.0_R8
      b2sihs_str=0.0_R8
      b2npmo_smaf=0.0_R8
      b2npmo_smag=0.0_R8
      b2npmo_smav=0.0_R8
      smpr=0.0_R8
      smpt=0.0_R8
      smfr=0.0_R8
      smcf=0.0_R8
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
      deallocate(smpr)
      deallocate(smpt)
      deallocate(smfr)
      deallocate(smcf)
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
      deallocate(b2stel_sna_ion)
      deallocate(b2stel_sna_rec)
      deallocate(b2stel_smq_ion)
      deallocate(b2stel_smq_rec)
      deallocate(b2stel_she)
      deallocate(b2stel_shi_ion)
      deallocate(b2stel_shi_rec)
      deallocate(b2stcx_sna)
      deallocate(b2stcx_smq)
      deallocate(b2stcx_shi)
      deallocate(b2srsm_sna)
      deallocate(b2srsm_smo)
      deallocate(b2srsm_she)
      deallocate(b2srsm_shi)
      deallocate(b2srdt_sna)
      deallocate(b2srdt_smo)
      deallocate(b2srdt_she)
      deallocate(b2srdt_shi)
      deallocate(b2srdt_sch)
      deallocate(b2srst_sna)
      deallocate(b2srst_smo)
      deallocate(b2srst_she)
      deallocate(b2srst_shi)
      deallocate(b2srst_sch)
      deallocate(b2stbm_sna)
      deallocate(b2stbm_smo)
      deallocate(b2stbm_sne)
      deallocate(b2stbm_she)
      deallocate(b2stbm_shi)
      deallocate(b2stbm_sch)
      deallocate(b2sihs_divue)
      deallocate(b2sihs_divua)
      deallocate(b2sihs_exbe)
      deallocate(b2sihs_exba)
      deallocate(b2sihs_visa)
      deallocate(b2sihs_joule)
      deallocate(b2sihs_fraa)
      deallocate(b2sihs_str)
      deallocate(ext_sna)
      deallocate(ext_smo)
      deallocate(ext_sne)
      deallocate(ext_she)
      deallocate(ext_shi)
      deallocate(ext_sch)
      deallocate(schdt)                        !srv 16.07.10 {
      deallocate(snedt)
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
      if(b2fplasma_version.ge.'01.001.026')
     >  call cfrure_opt (nread, size(sne), sne, 'sne')
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
      if(b2fplasma_version.ge.'01.001.026') call cfrure_opt 
     >  (nread, size(b2stbr_sne), b2stbr_sne, 'b2stbr_sne')
      if(b2fplasma_version.gt.'01.001.010') then
        call cfrure (nread, size(b2stbc_sna), b2stbc_sna, 'b2stbc_sna')
        call cfrure (nread, size(b2stbc_smo), b2stbc_smo, 'b2stbc_smo')
        call cfrure (nread, size(b2stbc_she), b2stbc_she, 'b2stbc_she')
        call cfrure (nread, size(b2stbc_shi), b2stbc_shi, 'b2stbc_shi')
        call cfrure (nread, size(b2stbc_sch), b2stbc_sch, 'b2stbc_sch')
      endif
      if(b2fplasma_version.ge.'03.000.006') then !djm 14/10/16
        call cfrure (nread, size(b2stel_sna_ion),
     .                b2stel_sna_ion, 'b2stel_sna_ion')
        call cfrure (nread, size(b2stel_sna_rec),
     .                b2stel_sna_rec, 'b2stel_sna_rec')
        call cfrure (nread, size(b2stel_smq_ion),
     .                b2stel_smq_ion, 'b2stel_smq_ion')
        call cfrure (nread, size(b2stel_smq_rec),
     .                b2stel_smq_rec, 'b2stel_smq_rec')
        call cfrure (nread, size(b2stel_she),
     .                b2stel_she, 'b2stel_she')
        call cfrure (nread, size(b2stel_shi_ion),
     .                b2stel_shi_ion, 'b2stel_shi_ion')
        call cfrure (nread, size(b2stel_shi_rec),
     .                b2stel_shi_rec, 'b2stel_shi_rec')
        call cfrure (nread, size(b2stcx_sna), b2stcx_sna, 'b2stcx_sna')
        call cfrure (nread, size(b2stcx_smq), b2stcx_smq, 'b2stcx_smq')
        call cfrure (nread, size(b2stcx_shi), b2stcx_shi, 'b2stcx_shi')
        call cfrure (nread, size(b2srsm_sna), b2srsm_sna, 'b2srsm_sna')
        call cfrure (nread, size(b2srsm_smo), b2srsm_smo, 'b2srsm_smo')
        call cfrure (nread, size(b2srsm_she), b2srsm_she, 'b2srsm_she')
        call cfrure (nread, size(b2srsm_shi), b2srsm_shi, 'b2srsm_shi')
        call cfrure (nread, size(b2srdt_sna), b2srdt_sna, 'b2srdt_sna')
        call cfrure (nread, size(b2srdt_smo), b2srdt_smo, 'b2srdt_smo')
        call cfrure (nread, size(b2srdt_she), b2srdt_she, 'b2srdt_she')
        call cfrure (nread, size(b2srdt_shi), b2srdt_shi, 'b2srdt_shi')
        call cfrure (nread, size(b2srdt_sch), b2srdt_sch, 'b2srdt_sch')
        call cfrure (nread, size(b2srst_sna), b2srst_sna, 'b2srst_sna')
        call cfrure (nread, size(b2srst_smo), b2srst_smo, 'b2srst_smo')
        call cfrure (nread, size(b2srst_she), b2srst_she, 'b2srst_she')
        call cfrure (nread, size(b2srst_shi), b2srst_shi, 'b2srst_shi')
        call cfrure (nread, size(b2srst_sch), b2srst_sch, 'b2srst_sch')
      endif
      if(b2fplasma_version.ge.'01.001.026') call cfrure_opt 
     >  (nread, size(b2stbc_sne), b2stbc_sne, 'b2stbc_sne')
      if(b2fplasma_version.ge.'01.001.027') then
        call cfrure_opt 
     >   (nread, size(b2stbm_sna), b2stbm_sna, 'b2stbm_sna')
        call cfrure_opt 
     >   (nread, size(b2stbm_smo), b2stbm_smo, 'b2stbm_smo')
        call cfrure_opt 
     >   (nread, size(b2stbm_she), b2stbm_she, 'b2stbm_she')
        call cfrure_opt 
     >   (nread, size(b2stbm_shi), b2stbm_shi, 'b2stbm_shi')
        call cfrure_opt 
     >   (nread, size(b2stbm_sch), b2stbm_sch, 'b2stbm_sch')
        call cfrure_opt 
     >   (nread, size(b2stbm_sne), b2stbm_sne, 'b2stbm_sne')
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
      if(b2fplasma_version.gt.'01.001.030') then
         call cfrure_opt (nread, size(smpr),  smpr, 'smpr')
         call cfrure_opt (nread, size(smpt),  smpt, 'smpt')
         call cfrure_opt (nread, size(smfr),  smfr, 'smfr')
      endif
      if(b2fplasma_version.gt.'03.000.000') then
         call cfrure (nread, size(smcf),  smcf, 'smcf')
      endif
      if(b2fplasma_version.gt.'01.001.012') then
        call cfrure (nread, size(ext_sna), ext_sna, 'ext_sna')
        call cfrure (nread, size(ext_smo), ext_smo, 'ext_smo')
        call cfrure (nread, size(ext_she), ext_she, 'ext_she')
        call cfrure (nread, size(ext_shi), ext_shi, 'ext_shi')
        call cfrure (nread, size(ext_sch), ext_sch, 'ext_sch')
      endif
      if(b2fplasma_version.ge.'01.001.026')
     >  call cfrure_opt (nread, size(ext_sne), ext_sne, 'ext_sne')

      return
      end subroutine read_b2mod_sources

      subroutine write_b2mod_sources(nwrite)
      implicit none
      integer nwrite
      character*3 source_name

      if (newversion.ge.'01.001.026') then
        source_name = 'sne'
      else if (newversion.gt.'01.001.010') then
        source_name = 'sch'
      end if

      call cfwure (nwrite, size(sch), sch, 'sch')
      call cfwure (nwrite, size(she), she, 'she')
      call cfwure (nwrite, size(shi), shi, 'shi')
      call cfwure (nwrite, size(smo), smo, 'smo')
      call cfwure (nwrite, size(smq), smq, 'smq')
      call cfwure (nwrite, size(sna), sna, 'sna')
      if(newversion.ge.'01.001.026')
     >  call cfwure (nwrite, size(sne), sne, 'sne')
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
        write(*,*) newversion,'b2stbr_sna --- _'//source_name
        call cfwure (nwrite, size(b2stbr_sna), b2stbr_sna, 'b2stbr_sna')
        call cfwure (nwrite, size(b2stbr_smo), b2stbr_smo, 'b2stbr_smo')
        call cfwure (nwrite, size(b2stbr_she), b2stbr_she, 'b2stbr_she')
        call cfwure (nwrite, size(b2stbr_shi), b2stbr_shi, 'b2stbr_shi')
        call cfwure (nwrite, size(b2stbr_sch), b2stbr_sch, 'b2stbr_sch')
      endif
      if(newversion.ge.'01.001.026')
     >  call cfwure (nwrite, size(b2stbr_sne), b2stbr_sne, 'b2stbr_sne')
      if(newversion.gt.'01.001.010') then
        write(*,*) newversion,'b2stbc_sna --- _'//source_name
        call cfwure (nwrite, size(b2stbc_sna), b2stbc_sna, 'b2stbc_sna')
        call cfwure (nwrite, size(b2stbc_smo), b2stbc_smo, 'b2stbc_smo')
        call cfwure (nwrite, size(b2stbc_she), b2stbc_she, 'b2stbc_she')
        call cfwure (nwrite, size(b2stbc_shi), b2stbc_shi, 'b2stbc_shi')
        call cfwure (nwrite, size(b2stbc_sch), b2stbc_sch, 'b2stbc_sch')
      endif
      if(newversion.ge.'03.000.006') then !djm 14/10/16
        call cfwure (nwrite, size(b2stel_sna_ion),
     .                b2stel_sna_ion, 'b2stel_sna_ion')
        call cfwure (nwrite, size(b2stel_sna_rec),
     .                b2stel_sna_rec, 'b2stel_sna_rec')
        call cfwure (nwrite, size(b2stel_smq_ion),
     .                b2stel_smq_ion, 'b2stel_smq_ion')
        call cfwure (nwrite, size(b2stel_smq_rec),
     .                b2stel_smq_rec, 'b2stel_smq_rec')
        call cfwure (nwrite, size(b2stel_she),
     .                b2stel_she, 'b2stel_she')
        call cfwure (nwrite, size(b2stel_shi_ion),
     .                b2stel_shi_ion, 'b2stel_shi_ion')
        call cfwure (nwrite, size(b2stel_shi_rec),
     .                b2stel_shi_rec, 'b2stel_shi_rec')
        call cfwure (nwrite, size(b2stcx_sna), b2stcx_sna, 'b2stcx_sna')
        call cfwure (nwrite, size(b2stcx_smq), b2stcx_smq, 'b2stcx_smq')
        call cfwure (nwrite, size(b2stcx_shi), b2stcx_shi, 'b2stcx_shi')
        call cfwure (nwrite, size(b2srsm_sna), b2srsm_sna, 'b2srsm_sna')
        call cfwure (nwrite, size(b2srsm_smo), b2srsm_smo, 'b2srsm_smo')
        call cfwure (nwrite, size(b2srsm_she), b2srsm_she, 'b2srsm_she')
        call cfwure (nwrite, size(b2srsm_shi), b2srsm_shi, 'b2srsm_shi')
        call cfwure (nwrite, size(b2srdt_sna), b2srdt_sna, 'b2srdt_sna')
        call cfwure (nwrite, size(b2srdt_smo), b2srdt_smo, 'b2srdt_smo')
        call cfwure (nwrite, size(b2srdt_she), b2srdt_she, 'b2srdt_she')
        call cfwure (nwrite, size(b2srdt_shi), b2srdt_shi, 'b2srdt_shi')
        call cfwure (nwrite, size(b2srdt_sch), b2srdt_sch, 'b2srdt_sch')
        call cfwure (nwrite, size(b2srst_sna), b2srst_sna, 'b2srst_sna')
        call cfwure (nwrite, size(b2srst_smo), b2srst_smo, 'b2srst_smo')
        call cfwure (nwrite, size(b2srst_she), b2srst_she, 'b2srst_she')
        call cfwure (nwrite, size(b2srst_shi), b2srst_shi, 'b2srst_shi')
        call cfwure (nwrite, size(b2srst_sch), b2srst_sch, 'b2srst_sch')
      endif
      if(newversion.ge.'01.001.026')
     >  call cfwure (nwrite, size(b2stbc_sne), b2stbc_sne, 'b2stbc_sne')
      if(newversion.ge.'01.001.027') then
        write(*,*) newversion,'b2stbm_sna --- _'//source_name
        call cfwure (nwrite, size(b2stbm_sna), b2stbm_sna, 'b2stbm_sna')
        call cfwure (nwrite, size(b2stbm_smo), b2stbm_smo, 'b2stbm_smo')
        call cfwure (nwrite, size(b2stbm_she), b2stbm_she, 'b2stbm_she')
        call cfwure (nwrite, size(b2stbm_shi), b2stbm_shi, 'b2stbm_shi')
        call cfwure (nwrite, size(b2stbm_sch), b2stbm_sch, 'b2stbm_sch')
        call cfwure (nwrite, size(b2stbm_sne), b2stbm_sne, 'b2stbm_sne')
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
        call cfwure (nwrite, size(b2npmo_smaf), b2npmo_smaf,
     1   'b2npmo_smaf')
        call cfwure (nwrite, size(b2npmo_smag), b2npmo_smag,
     1   'b2npmo_smag')
        call cfwure (nwrite, size(b2npmo_smav), b2npmo_smav,
     1   'b2npmo_smav')
      endif
      if(newversion.gt.'01.001.030') then
         call cfwure (nwrite, size(smpr),  smpr, 'smpr')
         call cfwure (nwrite, size(smpt),  smpt, 'smpt')
         call cfwure (nwrite, size(smfr),  smfr, 'smfr')
      endif
      if(newversion.gt.'03.000.000') then
         call cfwure (nwrite, size(smcf),  smcf, 'smcf')
      endif
      if(newversion.gt.'01.001.012') then
        call cfwure (nwrite, size(ext_sna), ext_sna, 'ext_sna')
        call cfwure (nwrite, size(ext_smo), ext_smo, 'ext_smo')
        call cfwure (nwrite, size(ext_she), ext_she, 'ext_she')
        call cfwure (nwrite, size(ext_shi), ext_shi, 'ext_shi')
        call cfwure (nwrite, size(ext_sch), ext_sch, 'ext_sch')
      endif
      if(newversion.ge.'01.001.026')
     >  call cfwure (nwrite, size(ext_sne), ext_sne, 'ext_sne')

      return
      end subroutine write_b2mod_sources

      end module b2mod_sources
