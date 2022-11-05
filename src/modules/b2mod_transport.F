      module b2mod_transport
      use b2mod_types
      use b2mod_version
      implicit none
      real (kind=R8), allocatable, save ::
     *  calf(:,:,:,:),
     *  cdna(:,:,:,:,:),
     *  cdpa(:,:,:,:,:),
     *  cddi(:,:,:,:,:),                                                !srv 13.10.06
     *  cdnahz(:,:,:,:,:),                                              !srv 09.01.01
     *  cdpahz(:,:,:,:,:),                                              !srv 09.01.01
     *  cdkt(:,:,:,:),
     *  cdzt(:,:,:,:),
     *  ceqp(:,:),
     *  chce(:,:,:,:),
     *  chci(:,:,:,:),
     *  chve(:,:,:,:),
     *  chvemx(:,:),
     *  chvi(:,:,:,:),
     *  chvimx(:,:),
     *  csig(:,:,:,:),
     *  cvla(:,:,:,:,:),
     *  cvsa(:,:,:,:,:),
     *  cthe(:,:,:),
     *  cthi(:,:,:),
     *  csigin(:,:,:,:,:,:),                                            !xpb
     *  cvlahz(:,:,:,:,:),                                              !srv 09.01.01
     *  cvmahz(:,:,:,:,:),                                              !srv 15.12.05
     *  cvsahz(:,:,:,:,:),                                              !srv 09.01.01
     *  cvsa_cl(:,:,:,:,:),                                             !srv 04.10.99
     *  cvsahz_cl(:,:,:,:,:),                                           !srv 09.01.01
     *  csig_cl(:,:,:,:),
     *  calf_cl(:,:,:,:),
     *  fllime(:,:),
     *  fllimi(:,:),
     *  fllim0fna(:,:,:,:,:),
     *  fllim0fhi(:,:,:,:,:),
     *  fllimvisc(:,:,:),
     *  fllim_ki(:,:),                                                  !srv 13.01.17 {
     *  fllim_ke(:,:),
     *  fllim_al(:,:),
     *  fllim_al_c(:,:),
     *  f_luc_ke(:,:),
     *  f_luc_ki(:,:),
     *  f_luc_et(:,:),
     *  f_luc_sg(:,:),
     *  f_luc_al(:,:),
     *  alfx_c(:,:),
     *  sigx_c(:,:),                                                    !iyv 19.11.13 08.10.14 !srv 13.01.17 }
     *  kt_neo(:,:)                                                     !iys 02.12.13 08.10.14

      contains

      subroutine alloc_b2mod_transport(nxd,nyd,nsd)
      implicit none
      integer nxd,nyd,nsd

      allocate(calf(-1:nxd,-1:nyd,0:1,0:1))
      allocate(cdna(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))
      allocate(cdpa(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))
      allocate(cddi(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))                     !srv 13.10.06
      allocate(cdnahz(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))                   !srv 09.01.01
      allocate(cdpahz(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))                   !srv 09.01.01
      allocate(cdkt(-1:nxd,-1:nyd,0:1,0:1))
      cdkt = 0.0_R8
      allocate(cdzt(-1:nxd,-1:nyd,0:1,0:1))
      cdzt = 0.0_R8
      allocate(ceqp(-1:nxd,-1:nyd))
      allocate(chce(-1:nxd,-1:nyd,0:1,0:1))
      allocate(chci(-1:nxd,-1:nyd,0:1,0:1))
      allocate(chve(-1:nxd,-1:nyd,0:1,0:1))
      allocate(chvemx(-1:nxd,-1:nyd))
      allocate(chvi(-1:nxd,-1:nyd,0:1,0:1))
      allocate(chvimx(-1:nxd,-1:nyd))
      allocate(csig(-1:nxd,-1:nyd,0:1,0:1))
      allocate(cvla(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))
      allocate(cvsa(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))
      allocate(cthe(-1:nxd,-1:nyd,0:nsd-1))
      allocate(cthi(-1:nxd,-1:nyd,0:nsd-1))
      allocate(csigin(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1,0:nsd-1))           !xpb
      allocate(cvlahz(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))                   !srv 09.01.01
      allocate(cvmahz(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))                   !srv 15.12.05
      allocate(cvsahz(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))                   !srv 09.01.01
      allocate(cvsa_cl(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))                  !srv 04.10.99
      allocate(cvsahz_cl(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))                !srv 09.01.01
      allocate(csig_cl(-1:nxd,-1:nyd,0:1,0:1))
      allocate(calf_cl(-1:nxd,-1:nyd,0:1,0:1))
      allocate(fllime(-1:nxd,-1:nyd))
      fllime=0.0_R8
      allocate(fllimi(-1:nxd,-1:nyd))
      fllimi=0.0_R8
      allocate(fllim0fna(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))
      fllim0fna=0.0_R8
      allocate(fllim0fhi(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))
      fllim0fhi=0.0_R8
      allocate(fllimvisc(-1:nxd,-1:nyd,0:nsd-1))
      fllimvisc=0.0_R8
      allocate(fllim_ki(-1:nxd,-1:nyd))                                 !srv 13.01.17 {
      allocate(fllim_ke(-1:nxd,-1:nyd))
      allocate(fllim_al(-1:nxd,-1:nyd))
      allocate(fllim_al_c(-1:nxd,-1:nyd))
      allocate(f_luc_ke(-1:nxd,-1:nyd))
      allocate(f_luc_ki(-1:nxd,-1:nyd))
      allocate(f_luc_et(-1:nxd,-1:nyd))
      allocate(f_luc_sg(-1:nxd,-1:nyd))
      allocate(f_luc_al(-1:nxd,-1:nyd))
      allocate(sigx_c(-1:nxd,-1:nyd))                                   !iys 19.11.13
      allocate(alfx_c(-1:nxd,-1:nyd))                                   !iys 19.11.13 08.10.14 ! srv 13.01.17 }
      allocate(kt_neo(-1:nxd,-1:nyd))                                   !iys 02.12.13 08.10.14
      kt_neo = 0.0_R8

      return
      end subroutine alloc_b2mod_transport

      subroutine dealloc_b2mod_transport
      implicit none

      if (.not.allocated(calf)) return
      deallocate(calf)
      deallocate(cdna)
      deallocate(cdpa)
      deallocate(cddi)
      deallocate(cdnahz)
      deallocate(cdpahz)
      deallocate(ceqp)
      deallocate(chce)
      deallocate(chci)
      deallocate(chve)
      deallocate(chvemx)
      deallocate(chvi)
      deallocate(chvimx)
      deallocate(csig)
      deallocate(cvla)
      deallocate(cvsa)
      deallocate(cthe)
      deallocate(cthi)
      deallocate(csigin)
      deallocate(cvlahz)
      deallocate(cvmahz)
      deallocate(cvsahz)
      deallocate(cvsa_cl)
      deallocate(cvsahz_cl)
      deallocate(csig_cl)
      deallocate(calf_cl)
      deallocate(fllime)
      deallocate(fllimi)
      deallocate(fllim0fna)
      deallocate(fllim0fhi)
      deallocate(fllimvisc)
      deallocate(fllim_ki)                                 !srv 13.01.17 {
      deallocate(fllim_ke)
      deallocate(fllim_al)
      deallocate(fllim_al_c)
      deallocate(f_luc_ke)
      deallocate(f_luc_ki)
      deallocate(f_luc_et)
      deallocate(f_luc_sg)
      deallocate(f_luc_al)
      deallocate(sigx_c)                                   !iys 19.11.13
      deallocate(alfx_c)                                   !iys 19.11.13 08.10.14 !srv 13.01.17 }
      deallocate(kt_neo)
      deallocate(cdkt)
      deallocate(cdzt)

      return
      end subroutine dealloc_b2mod_transport

      subroutine read_b2mod_transport(nxd,nyd,nsd,nread)
      implicit none
      integer        :: nxd,nyd,nsd, nread
      real (kind=R8) :: cloc(-1:nxd,-1:nyd,0:1), 
     &  clocs(-1:nxd,-1:nyd,0:1,0:nsd-1),
     &  clocss(-1:nxd,-1:nyd,0:1,0:nsd-1,0:nsd-1)
      

      if (b2fplasma_version.ge.'03.001.000') then
        call cfrure (nread, size(calf), calf, 'calf')
        call cfrure (nread, size(cdna), cdna, 'cdna')
        call cfrure (nread, size(cdpa), cdpa, 'cdpa')
        call cfrure (nread, size(ceqp), ceqp, 'ceqp')
        call cfrure (nread, size(chce), chce, 'chce')
        call cfrure (nread, size(chci), chci, 'chci')
        call cfrure (nread, size(chve), chve, 'chve')
        call cfrure (nread, size(chvemx), chvemx, 'chvemx')
        call cfrure (nread, size(chvi), chvi, 'chvi')
        call cfrure (nread, size(chvimx), chvimx, 'chvimx')
        call cfrure (nread, size(csig), csig, 'csig')
        call cfrure (nread, size(cvla), cvla, 'cvla')
        call cfrure (nread, size(cvsa), cvsa, 'cvsa')
      else
        call cfrure (nread, size(cloc), cloc, 'calf')
        calf = 0.0_R8
        calf(:,:,0,0) = cloc(:,:,0)
        calf(:,:,1,1) = cloc(:,:,1)
        call cfrure (nread, size(clocs), clocs, 'cdna')
        cdna = 0.0_R8
        cdna(:,:,0,0,:) = clocs(:,:,0,:)
        cdna(:,:,1,1,:) = clocs(:,:,1,:)
        call cfrure (nread, size(clocs), clocs, 'cdpa')
        cdpa = 0.0_R8
        cdpa(:,:,0,0,:) = clocs(:,:,0,:)
        cdpa(:,:,1,1,:) = clocs(:,:,1,:)
        call cfrure (nread, size(ceqp), ceqp, 'ceqp')
        call cfrure (nread, size(cloc), cloc, 'chce')
        chce = 0.0_R8
        chce(:,:,0,0) = cloc(:,:,0)
        chce(:,:,1,1) = cloc(:,:,1)
        call cfrure (nread, size(cloc), cloc, 'chci')
        chci = 0.0_R8
        chci(:,:,0,0) = cloc(:,:,0)
        chci(:,:,1,1) = cloc(:,:,1)
        call cfrure (nread, size(cloc), cloc, 'chve')
        chve = 0.0_R8
        chve(:,:,0,0) = cloc(:,:,0)
        chve(:,:,1,1) = cloc(:,:,1)
        call cfrure (nread, size(chvemx), chvemx, 'chvemx')
        call cfrure (nread, size(cloc), cloc, 'chvi')
        chvi = 0.0_R8
        chvi(:,:,0,0) = cloc(:,:,0)
        chvi(:,:,1,1) = cloc(:,:,1)
        call cfrure (nread, size(chvimx), chvimx, 'chvimx')
        call cfrure (nread, size(cloc), cloc, 'csig')
        csig = 0.0_R8
        csig(:,:,0,0) = cloc(:,:,0)
        csig(:,:,1,1) = cloc(:,:,1)
        call cfrure (nread, size(clocs), clocs, 'cvla')
        cvla = 0.0_R8
        cvla(:,:,0,0,:) = clocs(:,:,0,:)
        cvla(:,:,1,1,:) = clocs(:,:,1,:)
        call cfrure (nread, size(clocs), clocs, 'cvsa')
        cvsa = 0.0_R8
        cvsa(:,:,0,0,:) = clocs(:,:,0,:)
        cvsa(:,:,1,1,:) = clocs(:,:,1,:)
      endif
      call cfrure (nread, size(cthe), cthe, 'cthe')
      call cfrure (nread, size(cthi), cthi, 'cthi')
      if(b2fplasma_version.ge.'03.001.000') then
        call cfrure (nread, size(csigin), csigin, 'csigin')
        call cfrure (nread, size(cvsa_cl), cvsa_cl, 'cvsa_cl')
      elseif(b2fplasma_version.ge.'01.001.007') then
        call cfrure (nread, size(clocss), clocss, 'csigin')
        csigin = 0.0_R8
        csigin(:,:,0,0,:,:) = clocss(:,:,0,:,:)
        csigin(:,:,1,1,:,:) = clocss(:,:,1,:,:)
        call cfrure (nread, size(clocs), clocs, 'cvsa_cl')
        cvsa_cl = 0.0_R8
        cvsa_cl(:,:,0,0,:) = clocs(:,:,0,:)
        cvsa_cl(:,:,1,1,:) = clocs(:,:,1,:)
      else
        csigin=0.0_R8
        cvsa_cl=0.0_R8
      endif
      if(b2fplasma_version.ge.'01.001.021') then
        call cfrure (nread, size(fllime), fllime, 'fllime')
        call cfrure (nread, size(fllimi), fllimi, 'fllimi')
      else
        fllime=0.0_R8
        fllimi=0.0_R8
      endif
      if(b2fplasma_version.ge.'03.001.000') then
        call cfrure (nread, size(fllim0fna), fllim0fna, 'fllim0fna')
        call cfrure (nread, size(fllim0fhi), fllim0fhi, 'fllim0fhi')
      elseif(b2fplasma_version.ge.'01.001.017') then
        call cfrure (nread, size(clocs), clocs, 'fllim0fna')
        fllim0fna = 0.0_R8
        fllim0fna(:,:,0,0,:) = clocs(:,:,0,:)
        fllim0fna(:,:,1,1,:) = clocs(:,:,1,:)
        call cfrure (nread, size(clocs), clocs, 'fllim0fhi')
        fllim0fhi = 0.0_R8
        fllim0fhi(:,:,0,0,:) = clocs(:,:,0,:)
        fllim0fhi(:,:,1,1,:) = clocs(:,:,1,:)
      else
        fllim0fna=0.0_R8
        fllim0fhi=0.0_R8
      endif
      if(b2fplasma_version.ge.'01.001.020') then
        call cfrure (nread, size(fllimvisc), fllimvisc, 'fllimvisc')
      else
        fllimvisc=0.0_R8
      endif
      if(b2fplasma_version.ge.'03.000.006') then
        call cfrure (nread, size(f_luc_ke), f_luc_ke, 'f_luc_ke')
        call cfrure (nread, size(f_luc_ki), f_luc_ki, 'f_luc_ki')
        call cfrure (nread, size(f_luc_et), f_luc_et, 'f_luc_et')
        call cfrure (nread, size(f_luc_sg), f_luc_sg, 'f_luc_sg')
        call cfrure (nread, size(f_luc_al), f_luc_al, 'f_luc_al')
        call cfrure (nread, size(fllim_ke), fllim_ke, 'fllim_ke')
        call cfrure (nread, size(fllim_ki), fllim_ki, 'fllim_ki')
        call cfrure (nread, size(fllim_al), fllim_al, 'fllim_al')
        call cfrure (nread, size(fllim_al_c), fllim_al_c, 'fllim_al_c')
      else
        f_luc_ke=0.0_R8
        f_luc_ki=0.0_R8
        f_luc_et=0.0_R8
        f_luc_sg=0.0_R8
        f_luc_al=0.0_R8
        fllim_ke=0.0_R8
        fllim_ki=0.0_R8
        fllim_al=0.0_R8
        fllim_al_c=0.0_R8
      endif
      if(b2fplasma_version.ge.'03.002.000') then
        call cfrure (nread, size(cdkt), cdkt, 'cdkt')
      else
        cdkt=0.0_R8
      endif
      if(b2fplasma_version.ge.'03.002.000') then
        call cfrure (nread, size(cdzt), cdzt, 'cdzt')
      else
        cdzt=0.0_R8
      endif

      return
      end subroutine read_b2mod_transport

      subroutine write_b2mod_transport(nwrite)
      implicit none
      integer nwrite

      call cfwure (nwrite, size(calf), calf, 'calf')
      call cfwure (nwrite, size(cdna), cdna, 'cdna')
      call cfwure (nwrite, size(cdpa), cdpa, 'cdpa')
      call cfwure (nwrite, size(ceqp), ceqp, 'ceqp')
      call cfwure (nwrite, size(chce), chce, 'chce')
      call cfwure (nwrite, size(chci), chci, 'chci')
      call cfwure (nwrite, size(chve), chve, 'chve')
      call cfwure (nwrite, size(chvemx), chvemx, 'chvemx')
      call cfwure (nwrite, size(chvi), chvi, 'chvi')
      call cfwure (nwrite, size(chvimx), chvimx, 'chvimx')
      call cfwure (nwrite, size(csig), csig, 'csig')
      call cfwure (nwrite, size(cvla), cvla, 'cvla')
      call cfwure (nwrite, size(cvsa), cvsa, 'cvsa')
      call cfwure (nwrite, size(cthe), cthe, 'cthe')
      call cfwure (nwrite, size(cthi), cthi, 'cthi')
      if(newversion.ge.'01.001.007') then
        call cfwure (nwrite, size(csigin), csigin, 'csigin')
        call cfwure (nwrite, size(cvsa_cl), cvsa_cl, 'cvsa_cl')
      endif
      if(newversion.ge.'01.001.021') then
        call cfwure (nwrite, size(fllime), fllime, 'fllime')
        call cfwure (nwrite, size(fllimi), fllimi, 'fllimi')
      endif
      if(newversion.ge.'01.001.017') then
        call cfwure (nwrite, size(fllim0fna), fllim0fna, 'fllim0fna')
        call cfwure (nwrite, size(fllim0fhi), fllim0fhi, 'fllim0fhi')
      endif
      if(newversion.ge.'01.001.020') then
        call cfwure (nwrite, size(fllimvisc), fllimvisc, 'fllimvisc')
      endif
      if(newversion.ge.'03.000.006') then
        call cfwure (nwrite, size(f_luc_ke), f_luc_ke, 'f_luc_ke')
        call cfwure (nwrite, size(f_luc_ki), f_luc_ki, 'f_luc_ki')
        call cfwure (nwrite, size(f_luc_et), f_luc_et, 'f_luc_et')
        call cfwure (nwrite, size(f_luc_sg), f_luc_sg, 'f_luc_sg')
        call cfwure (nwrite, size(f_luc_al), f_luc_al, 'f_luc_al')
        call cfwure (nwrite, size(fllim_ke), fllim_ke, 'fllim_ke')
        call cfwure (nwrite, size(fllim_ki), fllim_ki, 'fllim_ki')
        call cfwure (nwrite, size(fllim_al), fllim_al, 'fllim_al')
        call cfwure (nwrite, size(fllim_al_c), fllim_al_c, 'fllim_al_c')
      endif
      if(newversion.ge.'03.002.000') then
        call cfwure (nwrite, size(cdkt), cdkt, 'cdkt')
      endif
      if(newversion.ge.'03.002.000') then
        call cfwure (nwrite, size(cdzt), cdzt, 'cdzt')
      endif

      return
      end subroutine write_b2mod_transport

      end module b2mod_transport
