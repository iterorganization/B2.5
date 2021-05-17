

      subroutine viscol(nsdecl, nzdecl, nsl, is, ic, ntau, m, t,
     +                  xi, ncl, x, den, nud, nue)
c--------------------------------------------------------------------
c     this subroutine calculates the energy-dependent collision
c     frequency and the energy-dependent energy scattering.
c
c     input  : nsdecl    the maximum number of species
c              nzdecl    the maximum number of charge states
c              nsl    the actual number of species
c              is     the species number for which the collision
c                     frequency is to be calculated
c              ic     the charge state of the species for which
c                     the collision frequency is to be calculated.
c              ntau   normalized collision frequency (n_i m_i /
c                     ntau_ij)  array(nsdecl+1,nsdecl+1)
c              m      array(nsdecl+1) the mass of the species in kg
c              t      array(nsdecl+1) the temperature of the species
c                     in kev
c              xi     array(nsdecl+1,nzdecl) the relative weight of every
c                     charge state
c              ncl    array(nsdecl+1) the number of charge states for 
c                     every species
c              x      normalized (to thermal) velocity for which
c                     the energy-dependent collision frequency is
c                     to be calculated
c              den    array(nsdecl+1,nzdecl) density of every component
c                     in 10^19 m^-3
c     output   nud    the pitch angle scattering frequency
c              nue    the energy scattering frequency
c
c     the routine calls the following function
c     erf   : calculates the error function (double precision)
c
c--------------------------------------------------------------------

      use b2mod_types

      implicit none

      integer ::
     *     nsdecl, nzdecl, is, ic, ncl(0:nsdecl), i, nsl
      real (kind=R8) ::
     *       ntau(0:nsdecl,0:nsdecl), m(0:nsdecl), t(0:nsdecl),
     &       xi(0:nsdecl,0:nzdecl-1), x, nud, nue,
     &       den(0:nsdecl,0:nzdecl-1)
      real (kind=R8) ::
     *       ph, g, pi, xab

      pi = 4.*atan(1.)

      nud = 0.
      nue = 0.

c     the loop over species
      do i = 0, nsl-1

c       calculate xab = vthb / vtha
        xab = sqrt(m(is)*t(i)/(m(i)*t(is)))

        ph = erf(x/xab)
        g  = (ph - 2*x*exp(-(x/xab)**2)/(xab*sqrt(pi)))
     +       /(2*(x/xab)**2)

        nud = nud + ntau(is,i) * ( ph - g )/ x**3
        nue = nue + ntau(is,i) * (-2.*ph/x**3 + 4.*(t(is)/t(i) +
     +        1./xab**2) *g/x)

      enddo

      nud = nud * 0.75 * sqrt(pi) * xi(is,ic) / (1.e19*den(is,ic)*
     +      m(is))
      nue = nue * 0.75 * sqrt(pi) * xi(is,ic) / (1.e19*den(is,ic)*
     +      m(is))

      return
      end subroutine viscol
