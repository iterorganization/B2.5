      module b2mod_fluid_molecules

      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_version
      use b2mod_constants
      use b2us_plasma
      use b2mod_subsys
      implicit none

      integer, save :: nms_ext
      integer, save :: background_was_read = 0

      real (kind=R8), allocatable, save ::
     .                             nm_ext(:,:), nmi_ext(:,:),
     .                             um_ext(:,:), umi_ext(:,:),
     .                             tm_ext(:,:), tmi_ext(:,:),
     .                             tau_mol(:), beta_eq(:), 
     .                             diff_coeff(:), coeff_vi(:),
     .                             diff_coeff_ene(:), 
     .                             uu_ext(:,:), vv_ext(:,:)
 
      real (kind=R8), parameter :: am_d2 = 4.0_R8 !! TODO generalize


      contains


      subroutine read_molecules_background_ext(nCv)
      implicit none

c input variables
      integer, intent(in) :: nCv
c local variables
      integer :: io, iun, idum(0:1), iex
      character :: check*10
#ifdef LEGACYCOMP
      integer newunit
      external newunit
#endif

      call subini ('read_molecules_background_ext')

#ifndef LEGACYCOMP
      open (newunit=iun,file='mol_bg',
     .    access='SEQUENTIAL',form='FORMATTED',status='OLD',iostat=iex)
#else
         iun = newunit()
      open (unit=iun,file='mol_bg',
     .    access='SEQUENTIAL',form='FORMATTED',status='OLD',iostat=iex)
#endif

      io = iex
      read (iun,'(A10)',iostat=io) check
*  read background of molecules and molecular ions provided by an external source
      if (io == 0) then
        backspace iun
        call cfruin (iun, 2, idum, 'nCv,nms')
        if (idum(0).ne.nCv) call xerrab('nCv not consistent')
        nms_ext = idum(1)
        if (.not.allocated(nm_ext)) then 
          allocate (nm_ext(idum(0),0:idum(1)-1)) !! molecules
          allocate (um_ext(idum(0),0:idum(1)-1))
          allocate (tm_ext(idum(0),0:idum(1)-1))
          allocate (nmi_ext(idum(0),0:idum(1)-1)) !! molecular ions
          allocate (umi_ext(idum(0),0:idum(1)-1))
          allocate (tmi_ext(idum(0),0:idum(1)-1))
          allocate (tau_mol(idum(0)))
          allocate (beta_eq(idum(0)))
          allocate (diff_coeff(idum(0)))
          allocate (coeff_vi(idum(0)))
          allocate (diff_coeff_ene(idum(0)))
          allocate (uu_ext(idum(0),0:idum(1)-1))
          allocate (vv_ext(idum(0),0:idum(1)-1))
        endif
        call cfrure (iun, idum(0)*idum(1), nm_ext, 'nm_ext')
        call cfrure (iun, idum(0)*idum(1), um_ext, 'um_ext')
        call cfrure (iun, idum(0)*idum(1), tm_ext, 'tm_ext')
        call cfrure (iun, idum(0)*idum(1), nmi_ext, 'nmi_ext')
        call cfrure (iun, idum(0)*idum(1), umi_ext, 'umi_ext')
        call cfrure (iun, idum(0)*idum(1), tmi_ext, 'tmi_ext')
        call cfrure (iun, idum(0)*idum(1), uu_ext, 'uu_ext')
        call cfrure (iun, idum(0)*idum(1), vv_ext, 'vv_ext')
        write (*,*) ' external molecule background read succesfully '
      else
        call xerrab('NO FILE mol_bg FOUND')
      end if

      call set_limits_mol_state(idum(0))
      background_was_read = 1

      close (iun)

      call subend ()
      return
      end subroutine read_molecules_background_ext

      subroutine set_limits_mol_state(nCv)
      use b2mod_types
      use b2mod_constants
      use b2mod_neutrals_namelist
     & , only : nm_min_ext
      implicit none
c input variables
      integer, intent(in) :: nCv
c local variables
      integer :: iCv
      real (kind=R8) :: nm_min, nm_max,
     .     um_max, tm_min, tm_max
      parameter ( nm_min=1.0e10_R8, nm_max=1.0e26_R8 )
      parameter ( um_max=1.0e5_R8 )   !! in absolute sense
      parameter ( tm_min=0.05_R8*ev, tm_max=1.0e4_R8*ev )


      do iCv = 1, nCv
        nm_ext(iCv,0) = min(nm_ext(iCv,0),nm_max)
        nm_ext(iCv,0) = max(nm_ext(iCv,0),nm_min_ext) !! make this one user defined for easy testing
        tm_ext(iCv,0) = min(tm_ext(iCv,0),tm_max)
        tm_ext(iCv,0) = max(tm_ext(iCv,0),tm_min)
        um_ext(iCv,0) =sign(
     &                  min(abs(um_ext(iCv,0)),um_max),
     &                          um_ext(iCv,0))
      enddo

      end subroutine set_limits_mol_state

      subroutine write_mol_state(iout,nCv)
      use b2mod_types
      use b2mod_constants
      implicit none
c input variables
      integer, intent(in) :: nCv, iout
c local variables
      integer :: nms, idum(0:1)

c to standard output
      call my_out_us(70,nCv,0,nm_ext(:,0),'nm_final')
      call my_out_us(70,nCv,0,um_ext(:,0),'um_final')
      call my_out_us(70,nCv,0,tm_ext(:,0),'tm_final')

c separate file to allow continuation
      idum(0) = nCv
      nms = 1 !! nms, TODO generalize
      idum(1) = nms
      
      call cfwuin (iout, 2, idum(0:1),'nCv,nms')

      call cfwure (iout, nCv*nms, nm_ext, 'nm_ext')
      call cfwure (iout, nCv*nms, um_ext, 'um_ext')
      call cfwure (iout, nCv*nms, tm_ext, 'tm_ext')
      call cfwure (iout, nCv*nms, nmi_ext, 'nmi_ext')
      call cfwure (iout, nCv*nms, umi_ext, 'umi_ext')
      call cfwure (iout, nCv*nms, tmi_ext, 'tmi_ext')
      call cfwure (iout, nCv*nms, uu_ext, 'uu_ext')
      call cfwure (iout, nCv*nms, vv_ext, 'vv_ext')

* ..return
      flush (iout)
      close (iout)
      return

      end subroutine write_mol_state

      end module b2mod_fluid_molecules
