      module b2mod_par_opt
      use b2mod_types

      integer :: npar_opt
      real (kind=R8), allocatable, save :: par_opt0(:), par_opt(:),
     &       par_opt_phys(:)
      logical, save :: optimization = .false.

      end module b2mod_par_opt
