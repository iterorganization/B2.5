c--------------------------------------------------------------------
c     start of the actual package
c--------------------------------------------------------------------

      subroutine neodv(nsdecl,nzdecl,nsl,ncl,zsp,ic, 
     +                 m,t,den,ds,rho,eparr,epsn,ddn,vvn,vtn)
c--------------------------------------------------------------------
c     subroutine that calculates the diffusion and pinch coefficients
c
c     input    : nsdecl   the maximum number of species
c                nzdecl   the maximum number of charge states
c                nsl    the actual number of species
c                ncl    array(nsl+1) the number of charge states per
c                       species
c                zsp   array(nsdecl+1,nzdecl) real. the charge of every
c                      component (normalized to e)
c                ic    the contribution for which the coefficients
c                      are calculated. 
c                      0 classical particle flux
c                      1 banana plateau contribution
c                      2 Pfirsch-Schlueter contribution
c                      3 both banana plateau and P-S
c                      4 all contributions
c                m     array(nsdecl+1) the mass of the species in kg
c                t     array(nsdecl+1) the temperature of the species in
c                      keV
c                den   array(nsdecl+1,nzdecl) the density of every component
c                      in units of 1.10^19 m^-3
c                ds    array(nsdecl+1,nzdecl,2) the thermodynamic forces 
c                      see the comments of the routine neoart to 
c                      find the exact definition
c                rho   the flux surface label
c                eparr loop voltage 
c                epsn  the desired accuracy
c     output     ddn   array(nsdecl+1,nzdecl) the diffusion coefficients 
c                      of every component in units of rho^2/sec
c                vvn   array(nsdecl+1,nzdecl) the total pinch coefficients of 
c                      every component in units of rho/sec
c                vtn   array(nsdecl+1,nzdecl) temperature dependent pinch  
c                      every component in units of rho/sec      
c
c     the routine calls the following subroutines
c     perr    : for error handling
c     neoart  : to calculate the fluxes
c--------------------------------------------------------------------
      use b2mod_types

      implicit none

      integer :: nsdecl
      integer :: nsl, ncl(0:nsdecl), nreg, nleg, nenergy, ncof, ic,
     +           i, j, k, l, mt, nzdecl
      real (kind=R8) :: m(0:nsdecl), t(0:nsdecl), 
     +       den(0:nsdecl,0:nzdecl-1),
     +       ds(0:nsdecl,0:nzdecl-1,2), dsl(0:nsdecl,0:nzdecl-1,2), 
     +       rho, epsn, 
     +       sigma(4), coeff(0:nsdecl,0:nzdecl-1,4), 
     +       ddn(0:nsdecl,0:nzdecl-1), vvn(0:nsdecl,0:nzdecl-1),
     +       zsp(0:nsdecl,0:nzdecl-1), vtn(0:nsdecl,0:nzdecl-1), 
     +       eparr, eparn
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

c     first initialize d and v to zero
      do  i = 0, nsl-1 
        do j = 0, ncl(i)-1
          ddn(i,j) = 0.
          vvn(i,j) = 0.
        enddo
      enddo

c     calculate the drift contribution
      do i=0,nsl-1
        do j=0,ncl(i)-1
            vtn(i,j) = 0.
        enddo
      enddo

c     set the epar to zero for the time being
      eparn = 0.
c     loop over the species 
      do i = 0, nsl-1
        do j = 0, ncl(i)-1

c       determine d
c       only normalized density gradient for i,j
        do k = 0, nsl-1
          do l = 0, ncl(k)-1
            do mt = 1, 2
              dsl(k,l,mt) = 0.
            enddo
          enddo
        enddo
        dsl(i,j,1) = 1.

        call neoart(nsl,ncl,nsdecl,nzdecl,zsp,m,t,den,dsl,rho,epsn,
     +              nreg,sigma,nleg,nenergy,ncof,
     +              neogeo,neofrc,ic,eparn,coeff)

c       the flux of component i,j appears in the diffusion
        ddn(i,j) = -coeff(i,j,1)/(den(i,j)*1e19)

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

        vvn(i,j) = coeff(i,j,1)/(den(i,j)*1e19)+ddn(i,j)*
     +            (ds(i,j,1)-ds(i,j,2))

        enddo
      enddo

      return 
      end subroutine neodv







