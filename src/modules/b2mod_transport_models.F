      module b2mod_transport_models
      use b2mod_types
      use b2mod_version
      logical, save :: transport_models_allocated=.false.
#ifdef ALLOCATE_AND_NAMELIST
      real (kind=R8), allocatable, save ::
     *  eta_hce_mult(:)
#else
#     include "DIMENSIONS.F"
      real (kind=R8), save :: eta_hce_mult(-1:DEF_NYD)
#endif

      namelist /transport_models/ eta_hce_mult

      contains

      subroutine alloc_b2mod_transport_models(nxd,nyd,nsd)
      implicit none
      integer nxd,nyd,nsd
      transport_models_allocated=.true.
#ifdef ALLOCATE_AND_NAMELIST
      allocate(eta_hce_mult(-1:nyd))
#else
      call xertst (nxd.le.DEF_NXD, 'Increase DEF_NXD in DIMENSIONS.F')
      call xertst (nyd.le.DEF_NYD, 'Increase DEF_NYD in DIMENSIONS.F')
      call xertst (nsd.le.DEF_NSD, 'Increase DEF_NSD in DIMENSIONS.F')
#endif
      end subroutine alloc_b2mod_transport_models

      subroutine dealloc_b2mod_transport_models
      implicit none
      if (.not.transport_models_allocated) return
#ifdef ALLOCATE_AND_NAMELIST
      deallocate(eta_hce_mult)
#endif
      transport_models_allocated=.false.
      end subroutine dealloc_b2mod_transport_models

      subroutine read_b2mod_transport_models(nread,edition,found)
      implicit none
      integer nread
      character*4 edition
      logical found
      logical :: file_ok
      character*260 filename
      external xertst

      call xertst (transport_models_allocated,
     .	'Error: transport_models not yet allocated')
      if (edition.eq.'    ') then
        filename = 'b2.transport_models_save.parameters'
      else
        filename = 'b2.transport_models_save.parameters.'//edition
      endif
      call find_file(filename,file_ok)
      if(file_ok) then
	if (edition.eq.'    ') open(nread,file=trim(filename))
	write(*,*) 'Reading namelist transport_models from ',
     1   trim(filename)
#ifdef USE_MPI
	read(nread,transport_models,err=999)
#else
        read(nread,transport_models)
#endif
	write(*,*) 'Read namelist transport_models from ',
     1   trim(filename)
        if (edition.eq.'    ') close(nread)
      else
	eta_hce_mult=1.0_R8
        write(6,*) 'No ',trim(filename)
      endif
      found = maxval(eta_hce_mult).ne.1.0_R8 .or.
     .        minval(eta_hce_mult).ne.1.0_R8
      return
 999  call xerrab ('Error while reading TRANSPORT_MODELS namelist!')
      end subroutine read_b2mod_transport_models

      subroutine write_b2mod_transport_models(nwrite,edition)
      implicit none
      integer nwrite
      character*4 edition
      character*260 filename

      if(.not.transport_models_allocated) return
      if (edition.eq.'    ') then
        filename = 'b2.transport_models_save.parameters'
      else
        filename = 'b2.transport_models_save.parameters.'//edition
      endif
      open(nwrite,file=trim(filename))
      write(nwrite,transport_models)
      close(nwrite)
      write(*,*) 'Saved ',trim(filename)
      return
      end subroutine write_b2mod_transport_models

      end module b2mod_transport_models
