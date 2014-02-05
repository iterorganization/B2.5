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
      contains

      subroutine alloc_b2mod_transport_models(nxd,nyd,nsd)
      integer nxd,nyd,nsd
      transport_models_allocated=.true.
#ifdef ALLOCATE_AND_NAMELIST
      allocate(eta_hce_mult(-1:nyd))
#else
      if(nyd.gt.DEF_NYD) stop 'Increase DEF_NYD in DIMENSIONS.F'
#endif
      end subroutine alloc_b2mod_transport_models

      subroutine read_b2mod_transport_models
      namelist /transport_models/ eta_hce_mult
      logical :: file_ok
      character*260 filename
      if(.not.transport_models_allocated) then
	write(*,*) 'Error: transport_models not yet allocated'
	stop
      endif
      filename='../' // 'b2.transport_models_save.parameters'
      inquire(file=filename,exist=file_ok)
      if(file_ok) then
	open(99,file=filename)
	write(*,*) 'Reading namelist transport_models from ',
     1	 trim(filename)
	read(99,transport_models)
	write(*,*) 'Read namelist transport_models from ',
     1	 trim(filename)
      else
	eta_hce_mult=1.0_R8
      endif
      end subroutine read_b2mod_transport_models

      subroutine write_b2mod_transport_models
      namelist /transport_models/ eta_hce_mult
      character*260 filename
      if(.not.transport_models_allocated) return
      filename='../' // 'b2.transport_models_save.parameters'
      open(99,file=filename)
      write(99,transport_models)
      write(*,*) 'Saved b2.transport_models_save.parameters'
      end subroutine write_b2mod_transport_models

      end module b2mod_transport_models
