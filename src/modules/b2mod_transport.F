      module b2mod_transport
      use b2mod_types
      use b2mod_version
      implicit none
      real (kind=R8), allocatable, save ::
     *  calf(:,:,:),
     *  cdna(:,:,:,:),
     *  cdpa(:,:,:,:),
     *  cddi(:,:,:,:),                                                  !srv 13.10.06
     *  cdde(:,:,:),                                                    !som 03.11.21
     *  cdnahz(:,:,:,:),                                                !srv 09.01.01
     *  cdpahz(:,:,:,:),                                                !srv 09.01.01
     *  ceqp(:,:),
     *  chce(:,:,:),
     *  chci(:,:,:),
     *  chve(:,:,:),
     *  chvemx(:,:),
     *  chvi(:,:,:),
     *  chvimx(:,:),
     *  csig(:,:,:),
     *  cvla(:,:,:,:),
     *  cvsa(:,:,:,:),
     *  cthe(:,:,:),
     *  cthi(:,:,:),
     *  csigin(:,:,:,:,:),                                              !xpb
     *  cvlahz(:,:,:,:),                                                !srv 09.01.01
     *  cvmahz(:,:,:,:),                                                !srv 15.12.05
     *  cvsahz(:,:,:,:),                                                !srv 09.01.01
     *  cvsa_cl(:,:,:,:),                                               !srv 04.10.99
     *  cvsa_drho(:,:,:,:),                                             !som 16.08.21
     *  cvsahz_cl(:,:,:,:),                                             !srv 09.01.01
     *  cvsahz_drho(:,:,:,:),                                           !som 18.08.21
     *  chci_a(:,:,:,:),
     *  csig_cl(:,:,:),
     *  calf_cl(:,:,:),
     *  csig_stoch(:,:,:),                                              !srv 24.08.21
     *  chce_stoch(:,:,:),                                              !srv 24.08.21
     *  fllime(:,:),
     *  fllimi(:,:),
     *  fllim0fna(:,:,:,:),
     *  fllim0fhi(:,:,:,:),
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
     *  alpha1(:,:),                                                    !som 05.10.21
     *  alfx_c(:,:),
     *  sigx_c(:,:),                                                    !iyv 19.11.13 08.10.14 !srv 13.01.17 }
     *  kt_neo(:,:)                                                     !iys 02.12.13 08.10.14

      contains

      subroutine alloc_b2mod_transport(nxd,nyd,nsd)
      implicit none
      integer, intent(in) :: nxd,nyd,nsd

      allocate(calf(-1:nxd,-1:nyd,0:1))
      allocate(cdna(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(cdpa(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(cddi(-1:nxd,-1:nyd,0:1,0:nsd-1))                         !srv 13.10.06
      allocate(cdde(-1:nxd,-1:nyd,0:1))                                 !som 03.11.21
      allocate(cdnahz(-1:nxd,-1:nyd,0:1,0:nsd-1))                       !srv 09.01.01
      allocate(cdpahz(-1:nxd,-1:nyd,0:1,0:nsd-1))                       !srv 09.01.01
      allocate(ceqp(-1:nxd,-1:nyd))
      allocate(chce(-1:nxd,-1:nyd,0:1))
      allocate(chci(-1:nxd,-1:nyd,0:1))
      allocate(chve(-1:nxd,-1:nyd,0:1))
      allocate(chvemx(-1:nxd,-1:nyd))
      allocate(chvi(-1:nxd,-1:nyd,0:1))
      allocate(chvimx(-1:nxd,-1:nyd))
      allocate(csig(-1:nxd,-1:nyd,0:1))
      allocate(cvla(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(cvsa(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(cthe(-1:nxd,-1:nyd,0:nsd-1))
      allocate(cthi(-1:nxd,-1:nyd,0:nsd-1))
      allocate(csigin(-1:nxd,-1:nyd,0:1,0:nsd-1,0:nsd-1))               !xpb
      csigin=0.0_R8
      allocate(cvlahz(-1:nxd,-1:nyd,0:1,0:nsd-1))                       !srv 09.01.01
      allocate(cvmahz(-1:nxd,-1:nyd,0:1,0:nsd-1))                       !srv 15.12.05
      allocate(cvsahz(-1:nxd,-1:nyd,0:1,0:nsd-1))                       !srv 09.01.01
      allocate(cvsa_cl(-1:nxd,-1:nyd,0:1,0:nsd-1))                      !srv 04.10.99
      allocate(cvsa_drho(-1:nxd,-1:nyd,0:1,0:nsd-1))                    !som 16.08.21
      cvsa_drho = 0.0_R8
      allocate(cvsahz_cl(-1:nxd,-1:nyd,0:1,0:nsd-1))                    !srv 09.01.01
      allocate(cvsahz_drho(-1:nxd,-1:nyd,0:1,0:nsd-1))                  !som 18.08.21
      cvsahz_drho = 0.0_R8      
      allocate(chci_a(-1:nxd,-1:nyd,0:1,0:nsd-1)) 
      chci_a = 0.0_R8
      allocate(csig_cl(-1:nxd,-1:nyd,0:1))
      allocate(calf_cl(-1:nxd,-1:nyd,0:1))
      allocate(csig_stoch(-1:nxd,-1:nyd,0:1))                           !srv 24.08.21
      allocate(chce_stoch(-1:nxd,-1:nyd,0:1))                           !srv 24.08.21
      csig_stoch = 0.0_R8
      chce_stoch = 0.0_R8
      allocate(fllime(-1:nxd,-1:nyd))
      fllime=1.0_R8
      allocate(fllimi(-1:nxd,-1:nyd))
      fllimi=1.0_R8
      allocate(fllim0fna(-1:nxd,-1:nyd,0:1,0:nsd-1))
      fllim0fna=1.0_R8
      allocate(fllim0fhi(-1:nxd,-1:nyd,0:1,0:nsd-1))
      fllim0fhi=1.0_R8
      allocate(fllimvisc(-1:nxd,-1:nyd,0:nsd-1))
      fllimvisc=1.0_R8
      allocate(fllim_ki(-1:nxd,-1:nyd))                                 !srv 13.01.17 {
      allocate(fllim_ke(-1:nxd,-1:nyd))
      allocate(fllim_al(-1:nxd,-1:nyd))
      allocate(fllim_al_c(-1:nxd,-1:nyd))
      allocate(f_luc_ke(-1:nxd,-1:nyd))
      allocate(f_luc_ki(-1:nxd,-1:nyd))
      allocate(f_luc_et(-1:nxd,-1:nyd))
      allocate(f_luc_sg(-1:nxd,-1:nyd))
      allocate(f_luc_al(-1:nxd,-1:nyd))
      allocate(alpha1(-1:nxd,-1:nyd))                                   !som 05.10.21
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
      deallocate(cdde)
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
      deallocate(cvsa_drho)                                      !som 16.08.21
      deallocate(cvsahz_cl)
      deallocate(cvsahz_drho)                                    !som 18.08.21
      deallocate(chci_a)
      deallocate(csig_cl)
      deallocate(calf_cl)
      deallocate(csig_stoch)                                     !srv 24.08.21
      deallocate(chce_stoch)                                     !srv 24.08.21
      deallocate(fllime)
      deallocate(fllimi)
      deallocate(fllim0fna)
      deallocate(fllim0fhi)
      deallocate(fllimvisc)
      deallocate(fllim_ki)                                       !srv 13.01.17 {
      deallocate(fllim_ke)
      deallocate(fllim_al)
      deallocate(fllim_al_c)
      deallocate(f_luc_ke)
      deallocate(f_luc_ki)
      deallocate(f_luc_et)
      deallocate(f_luc_sg)
      deallocate(f_luc_al)
      deallocate(alpha1)                                         !som 05.10.21
      deallocate(sigx_c)                                         !iys 19.11.13
      deallocate(alfx_c)                                         !iys 19.11.13 08.10.14 !srv 13.01.17 }
      deallocate(kt_neo)

      return
      end subroutine dealloc_b2mod_transport

      subroutine read_b2mod_transport(nread)
      implicit none
      integer, intent(in) :: nread
      external cfrure

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
      call cfrure (nread, size(cthe), cthe, 'cthe')
      call cfrure (nread, size(cthi), cthi, 'cthi')
      if(b2fplasma_version.ge.'01.001.007') then
        call cfrure (nread, size(csigin), csigin, 'csigin')
        call cfrure (nread, size(cvsa_cl), cvsa_cl, 'cvsa_cl')
      else
        csigin=0.0_R8
        cvsa_cl=0.0_R8
      endif
      if(b2fplasma_version.ge.'01.001.021') then
        call cfrure (nread, size(fllime), fllime, 'fllime')
        call cfrure (nread, size(fllimi), fllimi, 'fllimi')
      else
        fllime=1.0_R8
        fllimi=1.0_R8
      endif
      if(b2fplasma_version.ge.'01.001.017') then
        call cfrure (nread, size(fllim0fna), fllim0fna, 'fllim0fna')
        call cfrure (nread, size(fllim0fhi), fllim0fhi, 'fllim0fhi')
      else
        fllim0fna=1.0_R8
        fllim0fhi=1.0_R8
      endif
      if(b2fplasma_version.ge.'01.001.020') then
        call cfrure (nread, size(fllimvisc), fllimvisc, 'fllimvisc')
      else
        fllimvisc=1.0_R8
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
        f_luc_ke=1.0_R8
        f_luc_ki=1.0_R8
        f_luc_et=1.0_R8
        f_luc_sg=1.0_R8
        f_luc_al=1.0_R8
        fllim_ke=1.0_R8
        fllim_ki=1.0_R8
        fllim_al=1.0_R8
        fllim_al_c=1.0_R8
      endif

      return
      end subroutine read_b2mod_transport

      subroutine write_b2mod_transport(nwrite)
      implicit none
      integer, intent(in) :: nwrite
      external cfwure

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

      return
      end subroutine write_b2mod_transport

      end module b2mod_transport
