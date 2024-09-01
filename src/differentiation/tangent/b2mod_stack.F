      module b2mod_stack
      use b2mod_types

*     b2mod_stack contains two workspace stacks
      integer nstckw, nstcwi
      parameter (nstckw=10, nstcwi=10)
      real (kind=R8), save ::
     *   w(0:nstckw-1)
      integer, save :: wi(0:nstcwi-1)

      end module b2mod_stack
