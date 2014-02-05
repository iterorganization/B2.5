      module b2mod_b2cmwg
      use b2mod_types
*     (/b2cmwg/ contains some parameters to control graphical output)
      integer nmbpmx
      integer, save :: nmbpcc, nmbpdm
      parameter (nmbpmx=20)
      real (kind=R8), save ::
     *   cwaxr(0:nmbpmx-1), cwayr(0:nmbpmx-1), cwxmin(0:nmbpmx-1),
     &  cwxmax(0:nmbpmx-1), cwymin(0:nmbpmx-1), cwymax(0:nmbpmx-1)
      end module b2mod_b2cmwg
