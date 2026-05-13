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
     .                             n_mol(:,:), n_moli(:,:),
     .                             u_mol(:,:), u_moli(:,:),
     .                             t_mol(:,:), t_moli(:,:),
     .                             tau_mol(:), beta_eq(:), 
     .                             diff_coeff(:), coeff_vi(:),
     .                             diff_coeff_ene(:), 
     .                             u_tor_mol(:,:), u_dia_mol(:,:), !! toroidal, diagmagnetic
     .                             u_pol_mol(:,:), u_rad_mol(:,:)  !! poloidal, radial
 
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
        if (.not.allocated(n_mol)) then 
          allocate (n_mol(idum(0),0:idum(1)-1)) !! molecules
          allocate (u_mol(idum(0),0:idum(1)-1))
          allocate (t_mol(idum(0),0:idum(1)-1))
          allocate (n_moli(idum(0),0:idum(1)-1)) !! molecular ions
          allocate (u_moli(idum(0),0:idum(1)-1))
          allocate (t_moli(idum(0),0:idum(1)-1))
          allocate (tau_mol(idum(0)))
          allocate (beta_eq(idum(0)))
          allocate (diff_coeff(idum(0)))
          allocate (coeff_vi(idum(0)))
          allocate (diff_coeff_ene(idum(0)))
          allocate (u_tor_mol(idum(0),0:idum(1)-1))
          allocate (u_pol_mol(idum(0),0:idum(1)-1))
          allocate (u_dia_mol(idum(0),0:idum(1)-1))
          allocate (u_rad_mol(idum(0),0:idum(1)-1))
        endif
        call cfrure (iun, idum(0)*idum(1), n_mol, 'nm_ext')
        call cfrure (iun, idum(0)*idum(1), u_mol, 'um_ext')
        call cfrure (iun, idum(0)*idum(1), t_mol, 'tm_ext')
        call cfrure (iun, idum(0)*idum(1), n_moli, 'nmi_ext')
        call cfrure (iun, idum(0)*idum(1), u_moli, 'umi_ext')
        call cfrure (iun, idum(0)*idum(1), t_moli, 'tmi_ext')
        call cfrure (iun, idum(0)*idum(1), u_pol_mol, 'uu_ext')
        call cfrure (iun, idum(0)*idum(1), u_rad_mol, 'vv_ext')
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
        n_mol(iCv,0) = min(n_mol(iCv,0),nm_max)
        n_mol(iCv,0) = max(n_mol(iCv,0),nm_min_ext) !! make this one user defined for easy testing
        t_mol(iCv,0) = min(t_mol(iCv,0),tm_max)
        t_mol(iCv,0) = max(t_mol(iCv,0),tm_min)
        u_mol(iCv,0) =sign(
     &                  min(abs(u_mol(iCv,0)),um_max),
     &                          u_mol(iCv,0))
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
      call my_out_us(70,nCv,0,n_mol(:,0),'n_mol_final')
      call my_out_us(70,nCv,0,u_mol(:,0),'u_mol_final')
      call my_out_us(70,nCv,0,t_mol(:,0),'t_mol_final')

c separate file to allow continuation
      idum(0) = nCv
      nms = 1 !! nms, TODO generalize
      idum(1) = nms
      
      call cfwuin (iout, 2, idum(0:1),'nCv,nms')

      call cfwure (iout, nCv*nms, n_mol, 'nm_ext')
      call cfwure (iout, nCv*nms, u_mol, 'um_ext')
      call cfwure (iout, nCv*nms, t_mol, 'tm_ext')
      call cfwure (iout, nCv*nms, n_moli, 'nmi_ext')
      call cfwure (iout, nCv*nms, u_moli, 'umi_ext')
      call cfwure (iout, nCv*nms, t_moli, 'tmi_ext')
      call cfwure (iout, nCv*nms, u_pol_mol, 'uu_ext')
      call cfwure (iout, nCv*nms, u_rad_mol, 'vv_ext')

* ..return
      flush (iout)
      close (iout)
      return

      end subroutine write_mol_state

      end module b2mod_fluid_molecules
