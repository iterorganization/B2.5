c--------------------------------------------------------------------
c     start of the actual package
c--------------------------------------------------------------------

      subroutine neochi(nsdecl,nzdecl,nsl,ncl,zsp,ic,
     +                  m,t,den,ds,rho,eparr,epsn,kappan,chin,vwn)
c--------------------------------------------------------------------
c     subroutine that calculates the diffusion and pinch coefficients
c                      for the heat flux
c
c     input    : nsdecl   the maximum number of species (used for
c                         consistency check)
c                nzdecl   the maximum number of charge states (used
c                         for consistency check)
c                nsl    the actual number of species
c                ncl    array(nsdecl+1) the number of charge states per
c                       species
c                zsp   array(nsdecl+1,nzdecl) real. the charge of every
c                      component (normalized to e)
c                ic    the contribution for which the coefficients
c                      are calculated
c                      0 classical particle flux
c                      1 banana plateau contribution
c                      2 Pfirsch-Schlueter contribution
c                      3 both banana-plateau and P-S
c                      4 all contributions
c                m     array(nsdecl+1) the mass of the species in kg
c                t     array(nsdecl+1) the temperature of the species in
c                      keV
c                den   array(nsdecl+1,nzdecl) the density of every component
c                      in units of 10^19 m^-3
c                ds    array(nsdecl+1,nzdecl,2) the thermodynamic forces
c                      see the comments of the routine neoart to
c                      find the exact definition
c                rho   the flux surface label
c                eparr loop voltage
c                epsn   the desired accuracy
c     output     kappan array(nsdecl+1,nzdecl) the kappa coefficients
c                       of every component in units of rho^2/sec
c                chin   array(nsdecl+1,nzdecl) the chi coefficients
c                       of every component in units of rho^2/sec
c                vwn    array(nsdecl+1,nzdecl) the total pinch coefficients
c                       of every component in units of rho/sec
c
c     the routine calls the following subroutines
c     perr    : for error handling
c     neoart  : to calculate the fluxes
c
c     nnorm --- normalizes densities to 1e19 m^-3
c--------------------------------------------------------------------
      use b2mod_types

      implicit none

      real nnorm
      parameter(nnorm = 1e19)

      integer :: nsdecl
      integer :: nsl, ncl(0:nsdecl), nreg, nleg, nenergy, ncof, ic,
     +           i, j, k, l, mt, nzdecl
      real (kind=R8) :: m(0:nsdecl), t(0:nsdecl),
     +       den(0:nsdecl,0:nzdecl-1),
     +       ds(0:nsdecl,0:nzdecl-1,2), dsl(0:nsdecl,0:nzdecl-1,2),
     +       rho, epsn,
     +       sigma(4), coeff(0:nsdecl,0:nzdecl-1,4),
     +       kappan(0:nsdecl,0:nzdecl-1), zsp(0:nsdecl,0:nzdecl-1),
     +       eparr, eparn, chin(0:nsdecl,0:nzdecl-1),
     +       vwn(0:nsdecl,0:nzdecl-1)
      logical neofrc, neogeo

c     always use all the full expressions of the viscosity
      nreg = 0
c     use 3 legendre polynomials
      nleg = 3
c     use energy scattering in the collision operator
      nenergy = 1
c     use ion - electron collisions
      ncof = 1
c     use all couplings in the pfirsch schlueter regime
      sigma(1) = 1
      sigma(2) = 1
      sigma(3) = 1
      sigma(4) = 1
c     on the first call recalculate the friction/viscosity and
c     geometry parameters
      neofrc = .false.
      neogeo = .true.

c     first initialize kappan, chin and vwn to zero
      do i = 0, nsl-1
        do j = 0, ncl(i)-1
          kappan(i,j) = 0.
          chin(i,j) = 0.
          vwn(i,j) = 0.
        enddo
      enddo

c     set the epar to zero for the time being
      eparn = 0.
c     loop over the species
      do i = 0, nsl-1
        do j = 0, ncl(i)-1

c       determine kappan and chin
c       only normalized temperature gradient for i,j
        do k = 0, nsl-1
          do l = 0, ncl(k)-1
            do mt = 1, 2
              dsl(k,l,mt) = 0.
            enddo
          enddo
        enddo
        dsl(i,j,1) = 1.
        dsl(i,j,2) = 1.

        call neoart(nsl,ncl,nsdecl,nzdecl,zsp,m,t,den,dsl,rho,epsn,
     +              nreg,sigma,nleg,nenergy,ncof,
     +              neogeo,neofrc,ic,eparn,coeff)

c       the heat flux of component i,j appears in the diffusion
        kappan(i,j) = -coeff(i,j,2)
        chin(i,j) = kappan(i,j)/(den(i,j)*nnorm)

        enddo
      enddo

c     do not calculate a new friction and viscosity
      neogeo = .false.
      neofrc = .false.

c     calculate the pinch for the real thermodynamic forces
c       including the eparr
      call neoart(nsl,ncl,nsdecl,nzdecl,zsp,m,t,den,ds,rho,epsn,
     +            nreg,sigma,nleg,nenergy,ncof,
     +            neogeo,neofrc,ic,eparr,coeff)

c     loop over the species
      do i = 0, nsl-1
        do j = 0, ncl(i)-1

c       calculate pinch including the ware pinch

        vwn(i,j) = coeff(i,j,2) + kappan(i,j)*ds(i,j,2)

        enddo
      enddo

      return
      end subroutine neochi
