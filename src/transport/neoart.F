      subroutine neoart(nsp,ncs,nsd,nzd,zzs,msp,t,den,ds,rho,epsn,
     +                  nreg,sigma,nleg,nenergy,ncof,
     +                  neogeo,neofrc,ic,eparr,coeff)
c--------------------------------------------------------------------
c     subroutine this routine calculates the transport coefficients
c     of neoclassical theory.
c
c     input   nsp        : the number of species
c             ncs        : array(nsd+1) that gives the number of
c                          charge states per species
c             nsd        : actual (first) dimension of the zzs,den,
c                          and ds array
c             nzd        : actual (second) dimens. of array
c             zzs        : array(nsd+1,nzd) of real that contains
c                          the charge number of every component
c             msp        : array(nsd+1) that contains the mass of
c                          every species in kg
c             t          : array(nsdelc+1) that contains the temperature
c                          of the species in keV
c             den        : array(nsd+1,nzd) that contains the density
c                          of every component in units of 10^19 m^-3
c             ds         : array(nsd+1,nzd,2) that contains the first
c                          and second thermodynamic force (pressure
c                          gradient and temperature gradient)
c
c                                          d ln p_ij
c                              ds(i,j,1) = ---------
c                                            d rho
c
c                                          d ln t_i
c                              ds(i,j,2) = --------
c                                           d rho
c
c                          where rho is the flux label actually used
c                          the fluxes are also calculated in this
c                          coordinate. for this to be possible the
c                          geometry routine must supply dpsidr =
c                          d psi / d rho where psi is the poloidal
c                          flux.
c             rho        : the flux surface label
c             eps        : required accuracy
c             nreg       : parameter that forces a certain regime
c                          in the calculation of the viscosity
c                          0     the normal value. the banana and
c                                Pfirsch-Schlueter contributions
c                                are weighted
c                          1     the banana regime is forced.
c                          other the Pfirsch-Schlueter regime is
c                                forced.
c                          note : forcing the banana regime is not
c                                 the same thing as calculating the
c                                 banana plateau contribution. for
c                                 high collision freq. this contri-
c                                 bution usually decreases. when
c                                 nreg = 1 there is no such decr.
c             sigma(4)     the sigmas determine which of the coupling
c                          terms in the equation for the pfirsch
c                          schlueter regime are taken into acount.
c                          sigma(1)= 0,1  is the cross-coupling between
c                                     the heat equation and the equation
c                                     for the third laguerre harm.
c                          sigma(2)= 0.,1 is the cross-coupling between
c                                     the equation for the third harmonic
c                                     and the temperature perturb.
c                          sigma(3)= 0,1 is the energy exchange between
c                                     the different species
c                          sigma(4)= 0,1 is the coupling to the density
c                                     perturbation in the equation for the
c                                     third laguerre harmonic
c                          in the approximations of Hirshmann and Sigmar
c                          sigma(1)=0 sigma(2)=0 sigma(3)=1 sigma(4)=1
c                          that is the first two terms are neglected.
c                          the usual values of all these coefficients is 1.
c             nleg       : number of Legendre polynomials in the
c                          expansion (maximum and normal value is 3)
c             nenergy    : parameter that determines whether energy
c                          scattering is taken into account in the
c                          calculation of the viscosity
c                          0 not accounted for
c                          1 accounted for (normal value)
c             ncof         parameters that determines whether ion-
c                          electron collisions are taken into account
c                          0 not taken into account
c                          1 taken into account (normal value, 1;
c                          0 only to obtain some analytic results)
c             neogeo       logical if true then the geometry dep.
c                          parameters are recalculated
c             neofrc       logical if true then the friction and
c                          viscosity matrix are not newly calculated
c             ic         : the contribution for which the coefficients
c                          are calculated.
c                          0 the classical particle flux
c                          1 banana plateau contribution
c                          2 Pfirsch-Schlueter contribution
c                          3 both banana plateau, classical, and P-S
c                          4 all contributions
c             eparr      : the parallel electric field times the
c                          major radius, i.e. the loop voltage
c             coeff      : array(nsd+1,nzd,4) that contains the transp.
c                          coefficients of species nsm,ncm the last
c                          index gives
c                          1 particle flux   \gamma_ij
c                          2 energy flux     q_ij / t_i
c                          3 parallel flow   in units of current
c                                            <j.b><b>/<b**2>
c                          4 poloidal flow   <u.vt>/<b.vt>
c
c     the routine calls the following routines
c     perr    :  error handling
c     menn    :  calculate the friction coefficients
c     ps      :  calculate the pfirsch schlueter contribution
c     bp      :  calculate the banana plateau contribution
c--------------------------------------------------------------------

      use b2mod_types
      use b2mod_ranges
      use b2mod_constants
      use b2mod_neoclassical

      implicit none

      integer ::
     *       nsd,nzd,nreg,nleg,nenergy,ncof,ic
      real (kind=R8) ::
     *       ds(0:nsd,0:nzd-1,2),eparr
      logical neogeo, neofrc

      integer ::
     *       nsp,ncs(0:nsd),i,j,k,l
      real (kind=R8) ::
     *       msp(0:nsd),rho,den(0:nsd,0:nzd-1),t(0:nsd),
     &       ntau(0:nsd,0:nsd),coeff(0:nsd,0:nzd-1,4),
     &       zzs(0:nsd,0:nzd-1),eparn,epsn,
     &       xi(0:nsd,0:nzd-1),
     &       uai(0:nsd,0:nzd,3),sigma(4), dum,
     &       la(3,3,0:nsd), lab(3,3,0:nsd,0:nsd), mm(3,3),
     &       nn(3,3),ta,tb,ma,mb,
     &       coeffc(0:nsd,0:nzd-1,4)

      if  ((neofrc).and.(neogeo)) call perr(16)

c     if the friction and viscosity coefficients are newly cal-
c     culated then go through the following loop.
      if (.not.neofrc) then

      if (rho.eq.0.) call perr(2)

      if ((ic.eq.1).or.(ic.eq.3)) then
        if (nreg.ne.0) call perr(11)
        if (nleg.ne.3) call perr(12)
        if (nenergy.ne.1) call perr(13)
      endif
      if ((ic.eq.2).or.(ic.eq.3)) then
        dum = 1.
        do i = 1, 4
          dum = dum*sigma(i)
        enddo
        if (dum.ne.1.) call perr(15)
      endif
      if (ncof.ne.1) call perr(14)

      call colxi(nsd,nzd,nsp,ncs,zzs,den,t,msp,ntau,xi)

c     switch off ion - electron collisions?
      if (ncof.eq.0) then
c       yes switch off
        j = 0
        do i = 0, nsp-1
          if (zzs(i,0).lt.0.0_R8) j = i
        enddo
        if (j.ne.0) then
          do i = 0, nsp-1
            if (j.ne.i) ntau(i,j) = 0.0_R8
          enddo
        endif
      endif

c     calculate the friction coefficients
      do i = 0, nsp-1
        do k = 1, 3
          do l = 1, 3
            la(k,l,i) = 0.0_R8
          enddo
        enddo
        do j = 0, nsp-1
          ta = t(i)
          tb = t(j)
          ma = msp(i)
          mb = msp(j)
          call menn(ta,tb,ma,mb,mm,nn)
          do k = 1, 3
            do l = 1, 3
               lab(k,l,i,j) = ntau(i,j)*nn(k,l)
               la(k,l,i) = la(k,l,i) + ntau(i,j)*mm(k,l)
            enddo
          enddo
        enddo
      enddo

c     end of the if statement on neofrc
      endif

c     make the coefficients zero
      do i = 0, nsd
        do j = 0, nzd-1
          do k = 1,4
            coeff(i,j,k) = 0.0_R8
            coeffc(i,j,k) = 0.0_R8
          enddo
        enddo
      enddo

cc     normalize the eparr  (CHECK !!!)
      eparn = 1.0e-3_R8*dpsidr*r2i*eparr

c     calculate the classical transport contribution
      if ((ic.eq.0).or.(ic.eq.4)) then
        call class(nsd,nzd,nsp,ncs,t,zzs,xi,la,lab,ds,coeffc)
c       add the correct normalization factor, put result
c       in coeff and clear coeffc
        do i = 0, nsp-1
          do j = 0, ncs(i)-1
            do k = 1, 2
              coeff(i,j,k) = coeff(i,j,k)+ coeffc(i,j,k)*
     +          gclass/((ev/1000.0_R8)*zzs(i,j)*dpsidr**2)
              coeffc(i,j,k) = 0.0_R8
            enddo
          enddo
        enddo
      endif

c     calculate the Pfirsch-Schlueter contribution
      if ((ic.eq.2).or.(ic.eq.3).or.(ic.eq.4)) then
        call ps(nsd,nzd,nsp,ncs,neofrc,ntau,xi,la,lab,msp,t,den,
     +           zzs,sigma,ds,rnq,uai,coeffc)
c       add the correct normalization factor, put result in
c       coeff and clear coeffc
        do i = 0, nsp-1
          do j = 0, ncs(i)-1
            do k = 1, 2
              coeff(i,j,k) = coeff(i,j,k) - coeffc(i,j,k)*
     +          (1.- b2av*bi2a)*(rbt/dpsidr)**2/
     +          ((ev/1000.0_R8)*zzs(i,j)*b2av)
              coeffc(i,j,k) = 0.0_R8
            enddo
          enddo
        enddo
      endif

      if ((ic.eq.1).or.(ic.eq.3).or.(ic.eq.4)) then
        call bp(nsd,nzd,nsp,ncs,neofrc,t,msp,den,zzs,ntau,xi,la,
     +   lab,epsn,nenergy,nleg,nreg,bgradp,fc,fm,mmx,ds,eparn,
     +   coeffc)
c       add the correct normalizing factor. put result in
c       coeff and clear coeffc
        do i = 0, nsp-1
          do j = 0, ncs(i)-1
            do k = 1, 2
              coeff(i,j,k) = coeff(i,j,k) - coeffc(i,j,k)*
     +          (rbt/dpsidr)**2/((ev/1000.0_R8)*zzs(i,j)*b2av)
              coeffc(i,j,k) = 0.0_R8
            enddo
          enddo
        enddo
      endif

      do i = 0, nsp-1
        do j = 0, ncs(i)-1
          coeff(i,j,3) = coeffc(i,j,3)*(ev*1000.0_R8*1.0e19_R8)*
     +                   den(i,j)*zzs(i,j)*rbt*bav/(dpsidr*b2av)
          coeff(i,j,4) = 1.0e3_R8*coeffc(i,j,4)*rbt/(b2av*dpsidr)
        enddo
      enddo

      return
      end subroutine neoart



