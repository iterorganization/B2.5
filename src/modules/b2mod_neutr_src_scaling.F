      module b2mod_neutr_src_scaling

c  version : 24.07.2002 17:08

      use b2mod_types
      use b2mod_version
      real (kind=R8), allocatable, save ::
c      real sclstra(0:nsd,nstra),sclalph(0:nsd,nstra),
c     ,	 sclsrcn(0:nsd,nstra),sclbase(0:nsd,nstra),sclrtio(0:nsd,nstra)
     *  sclstra(:,:),
     *  sclalph(:,:),
     *  sclsrcn(:,:),
     *  sclbase(:,:),
     *  sclrtio(:,:),
     1	tflux(:),
     2	flux_scale(:)

      contains

      subroutine alloc_b2mod_neutr_src_scaling(nsd,nstra)
      implicit none
      integer nsd, nstra

      allocate(sclstra(0:nsd,nstra))
      allocate(sclalph(0:nsd,nstra))
      allocate(sclsrcn(0:nsd,nstra))
      allocate(sclbase(0:nsd,nstra))
      allocate(sclrtio(0:nsd,nstra))
      allocate(tflux(nstra))
      allocate(flux_scale(nstra))

      sclstra=1.0_R8
      sclalph=1.0_R8
      sclsrcn=0.0_R8
      sclbase=1.0_R8
      sclrtio=1.0_R8
      tflux=0.0_R8
      flux_scale=0.0_R8

      return
      end subroutine alloc_b2mod_neutr_src_scaling

      subroutine dealloc_b2mod_neutr_src_scaling
      implicit none

      if (.not.allocated(sclstra)) return
      deallocate(sclstra)
      deallocate(sclalph)
      deallocate(sclsrcn)
      deallocate(sclbase)
      deallocate(sclrtio)
      deallocate(tflux)
      deallocate(flux_scale)

      return
      end subroutine dealloc_b2mod_neutr_src_scaling

      subroutine read_b2mod_neutr_src_scaling(nread,ns,nstra)
      implicit none
#     include "DIMENSIONS.F"
      integer nread, ns, nstra, idum(0:9)

      if(b2fplasma_version.ge.'01.001.024') then
        call cfruin (nread, 1, idum, 'nstra')
        nstra=idum(0)
      else
        nstra=DEF_NSTRA
      endif
      call alloc_b2mod_neutr_src_scaling(ns,nstra)
      if(b2fplasma_version.ge.'01.001.023' .and. nstra.gt.0) then
        call cfrure (nread, size(sclstra), sclstra, 'sclstra')
        call cfrure (nread, size(sclrtio), sclrtio, 'sclrtio')
      endif

      return
      end subroutine read_b2mod_neutr_src_scaling

      subroutine write_b2mod_neutr_src_scaling(nwrite,nstra)
      implicit none
      integer nwrite, nstra, idum(0:9)

      if(newversion.ge.'01.001.024') then
        idum(0)=nstra
        call cfwuin (nwrite, 1, idum, 'nstra') 
      endif
      if(newversion.ge.'01.001.023' .and. nstra.gt.0) then
        write(*,*) newversion,'sclstra, sclrtio'
        call cfwure (nwrite, size(sclstra), sclstra, 'sclstra')
        call cfwure (nwrite, size(sclrtio), sclrtio, 'sclrtio')
      endif

      return
      end subroutine write_b2mod_neutr_src_scaling

      end module b2mod_neutr_src_scaling
