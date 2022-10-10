      module b2mod_par_opt
      use b2mod_types

      integer, save :: npar_opt = 0
      real (kind=R8), allocatable, save :: par_opt0(:), par_opt(:),
     &       par_opt_phys(:)

      end module b2mod_par_opt
