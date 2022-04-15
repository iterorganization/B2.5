      module b2mod_anomalous_transport
      use b2mod_types
      use b2mod_version
      implicit none
      real (kind=R8), allocatable, save ::
     *   sig0(:,:),
     *   hce0(:,:),
     *   alf0(:,:),
     *   hci0(:,:),
     *   hcib(:,:,:),
     *   chcb(:,:,:,:,:),
     *   dpa0(:,:,:),
     *   dna0(:,:,:),
     *   vsa0(:,:,:),
     *   vla0(:,:,:,:),
     *   csig_an(:,:,:,:),                                              !xpb
     *   calf_an(:,:,:,:),                                              !xpb
     *   cvsa_an(:,:,:,:,:),                                            !srv 04.10.99 11.06.08
     *   vma0(:,:,:,:)                                                  !srv 15.12.05 11.06.08

      contains

      subroutine dealloc_anomalous_transport
      implicit none

      deallocate(sig0)
      deallocate(hce0)
      deallocate(alf0)
      deallocate(hci0)
      deallocate(hcib)
      deallocate(chcb)
      deallocate(dpa0)
      deallocate(dna0)
      deallocate(vsa0)
      deallocate(vla0)
      deallocate(csig_an)
      deallocate(calf_an)
      deallocate(cvsa_an)                                               !srv 15.12.05 11.06.08
      deallocate(vma0)                                                  !srv 15.12.05 11.06.08

      return
      end subroutine dealloc_anomalous_transport

      subroutine alloc_b2mod_anomalous_transport(nxd,nyd,nsd)
      implicit none
      integer nxd,nyd,nsd                                               !xpb

      allocate(sig0(-1:nxd,-1:nyd))
      allocate(hce0(-1:nxd,-1:nyd))
      allocate(alf0(-1:nxd,-1:nyd))
      allocate(hci0(-1:nxd,-1:nyd))
      allocate(hcib(-1:nxd,-1:nyd,0:nsd-1))
      allocate(chcb(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))
      allocate(dpa0(-1:nxd,-1:nyd,0:nsd-1))
      allocate(dna0(-1:nxd,-1:nyd,0:nsd-1))
      allocate(vsa0(-1:nxd,-1:nyd,0:nsd-1))
      allocate(vla0(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(csig_an(-1:nxd,-1:nyd,0:1,0:1))                          !xpb
      allocate(calf_an(-1:nxd,-1:nyd,0:1,0:1))                          !xpb
      allocate(cvsa_an(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))                  !srv 04.10.99 11.06.08
      allocate(vma0(-1:nxd,-1:nyd,0:1,0:nsd-1))                         !srv 15.12.05 11.06.08

      return
      end subroutine alloc_b2mod_anomalous_transport

      subroutine read_b2mod_anomalous_transport(nxd,nyd,nread)
      implicit none
      integer :: nxd,nyd,nread
      real (kind=R8) :: cloc(-1:nxd,-1:nyd,0:1)
      external cfrure, cfrure_opt

      if(b2fplasma_version.gt.'01.001.001') then
        call cfrure (nread, size(sig0), sig0, 'sig0')
        call cfrure (nread, size(hce0), hce0, 'hce0')
        call cfrure (nread, size(alf0), alf0, 'alf0')
        call cfrure (nread, size(hci0), hci0, 'hci0')
        call cfrure (nread, size(hcib), hcib, 'hcib')
        call cfrure (nread, size(dpa0), dpa0, 'dpa0')
        call cfrure (nread, size(dna0), dna0, 'dna0')
        call cfrure (nread, size(vsa0), vsa0, 'vsa0')
        call cfrure (nread, size(vla0), vla0, 'vla0')
        if(b2fplasma_version.ge.'03.001.000') then
          call cfrure (nread, size(csig_an), csig_an, 'csig_an')
        else
          call cfrure (nread, size(cloc), cloc, 'csig_an')
          csig_an = 0.0_R8
          csig_an(:,:,0,0) = cloc(:,:,0)
          csig_an(:,:,1,1) = cloc(:,:,1)
        endif
      endif
      if(b2fplasma_version.ge.'03.001.000') then
        call cfrure_opt (nread, size(calf_an), calf_an, 'calf_an')
      elseif(b2fplasma_version.ge.'01.001.029') then
        call cfrure_opt (nread, size(cloc), cloc, 'calf_an')
        calf_an = 0.0_R8
        calf_an(:,:,0,0) = cloc(:,:,0)
        calf_an(:,:,1,1) = cloc(:,:,1)
      endif

      return
      end subroutine read_b2mod_anomalous_transport

      subroutine write_b2mod_anomalous_transport(nwrite)
      implicit none
      integer nwrite
      external cfwure

      if(newversion.gt.'01.001.001') then
        call cfwure (nwrite, size(sig0), sig0, 'sig0')
        call cfwure (nwrite, size(hce0), hce0, 'hce0')
        call cfwure (nwrite, size(alf0), alf0, 'alf0')
        call cfwure (nwrite, size(hci0), hci0, 'hci0')
        call cfwure (nwrite, size(hcib), hcib, 'hcib')
        call cfwure (nwrite, size(dpa0), dpa0, 'dpa0')
        call cfwure (nwrite, size(dna0), dna0, 'dna0')
        call cfwure (nwrite, size(vsa0), vsa0, 'vsa0')
        call cfwure (nwrite, size(vla0), vla0, 'vla0')
        call cfwure (nwrite, size(csig_an), csig_an, 'csig_an')
      endif
      if(newversion.ge.'01.001.029') then
        call cfwure (nwrite, size(calf_an), calf_an, 'calf_an')
      endif

      return
      end subroutine write_b2mod_anomalous_transport

      end module b2mod_anomalous_transport
