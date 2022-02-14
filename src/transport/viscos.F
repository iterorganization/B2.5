
      subroutine viscos(nsdecl, nzdecl, nsl, is, ic, nreg, 
     +              bgradp, fc, fm, mmx, nenergy, epsn, ntau, m, 
     +              t, xi, den, mu)
c--------------------------------------------------------------------
c     this routine calculates the viscosity coefficients using the 
c     expressions of K.C. Shaing, Phys. Plasmas 3 965 (1996)
c     extended to the 3 Laguerre polynomial expansion
c
c     input   nsdecl         : maximum number of species
c             nzdecl         : maximum number of charge states
c             nsl            : the actual number of species
c             is             : the species number for which the 
c                              viscosity is to be calculated
c             ic             : the charge state for which the visc-
c                              osity is to be calculated.
c             nreg           : integer that determines what regime
c                              is to be forced in the calculation 
c                              of the viscosity
c                              0 = weigted over all regimes
c                              1 = the banana regimes is forced
c                              other = the Pfirsch-Schlueter regime
c                              is forced. 
c             bgradp         : the flux surface quantity n in grad t
c                              (the inner product of the unit vector
c                              in the direction of the magnetic field
c                              and the gradient of the poloidal 
c                              angle)
c             fc             : the number of passing particles
c             fm             : the coefficients necessary to calculate
c                              the viscosity in the pfirsch schlueter
c                              regime.
c             mmx            : the number of coefficients fm
c             epsn           : the accuracy with which the numerical
c                              coeficients should be calculated. 
c                              note that the accuracy of the analytic
c                              equations implemented here is not 
c                              expected to be larger than 10% in 
c                              certain collisional regimes. when epsn = 0., 
c                              the value 0.001 will be assumed.
c             ntau           : array(nsdecl+1,nsdecl+1) the collision 
c                              frequency n_i m_i / ntau_ij
c             m              : array(nsdecl+1) the masses of the species 
c                              in kg
c             t              : array(nsdecl+1) temperature of the species 
c                              in kev
c             xi             : array(nsdecl+1,nzdecl) relative weigth of a 
c                              charge state
c             den            : array(nsdecl+1,nzdecl) density in 
c                              unit 10^19 m^-3
c     output  mu             : array(3,3) the viscosity coefficients
c
c     the routine calls the following subroutine
c     viscol    : routine that calculates the energy dependent 
c                 collision and energy scattering frequency
c     perr      : routine that does the error handling.
c--------------------------------------------------------------------
       
      use b2mod_types

      implicit none

      integer ::
     *     nenergy, nsdecl, nzdecl, is, ic, nsl, nreg
      real (kind=R8) ::
     *       epsn, mu(3,3), ntau(0:nsdecl,0:nsdecl), t(0:nsdecl), 
     &       m(0:nsdecl), xi(0:nsdecl,0:nzdecl-1), 
     &       den(0:nsdecl,0:nzdecl-1)
      integer ::
     *     mmx,i,j,nv,k,l
      logical nlerr
      real (kind=R8) ::
     *       bgradp, fc, ft, pi, x, dum, vth,
     +       nud, nue, nut, kb, kps, ommn, ntom, ktot, voorf,
     +       dx, weight, meanmuo
      real (kind=R8) ::
     *       fm(mmx), a(6), muo(3,3)
      real (kind=R8) ::
     *       xpoints(0:25600), cxi(0:25600)

c     calculate vth
      vth = sqrt(2.*1.6e-16*t(is)/m(is))

c     the constant pi, 2 pi
      pi = 4.*atan(1.)

c     the constants to approximate the nut*irm function
      a(1) = 2. / 5.
      a(2) = - 22. / 105.
      a(3) = 6. / 35.
      a(4) = - 34. / 231.
      a(5) = 166. / 1287.
      a(6) = - 82. / 715.


      ft = 1.- fc

      nv = 100
      do k = 1, 3
        do l = 1, 3
          muo(k,l) = 0.
        enddo
      enddo
 3000 continue

c       the loop over the velocity
        dx = 10./real(nv)
        do k = 1, 3
          do l = 1, 3
            mu(k,l) = 0.
          enddo
        enddo

*   ..integration via Simpson rule

        do i = 0, nv 

          x = i*dx
          ommn = vth * x * abs(bgradp)
        
c         now calculate the collision and energy scattering 
c         frequency. 

          if (i.eq.0) then

           voorf = 0.

          else

           call viscol(nsdecl, nzdecl, nsl, is, ic, ntau, m, t,
     +                 xi, x, den, nud, nue)
           nut = 3*nud + real(nenergy)*nue



           kb  = ft * nud / fc
           kps = 0.
           do j = 1, mmx
             ntom = nut / ( ommn * real(j) )
             if (ntom.lt.10) then
               dum  = -1.5*ntom**2 - 4.5*ntom**4 + 
     +              (0.25+(1.5+2.25*ntom**2)*ntom**2)*2*ntom *
     +              atan(1/ntom)
             else
               ntom = (1. / ntom)**2 
               dum = a(1)
               do k = 2, 6
                 dum = dum + a(k)*ntom**(k-1)
               enddo
             endif
             kps = kps + fm(j)*dum
           enddo
           kps = kps * 1.5 * (vth * x)**2/nut
           if (nreg.eq.0) then
            ktot = kb*kps/(kb+kps)
           else
             if (nreg.eq.1) then
               ktot = kb
             else 
               ktot = kps
             endif
           endif

           voorf = x**4 * exp(-x**2)*ktot*dx/3.

          endif

          if ((i.eq.0).or.(i.eq.nv)) then
             weight=1.
           else
             weight=2.*2.**mod(i,2)
          endif
          mu(1,1) = mu(1,1) + weight*voorf
          mu(1,2) = mu(1,2) + weight*voorf*(x**2-2.5)
          mu(2,2) = mu(2,2) + weight*voorf*(x**2-2.5)**2
          mu(1,3) = mu(1,3) + weight*voorf*(35./8.-3.5*x**2 
     +                        + 0.5*x**4)
          mu(2,3) = mu(2,3) + weight*voorf*(35./8.-3.5*x**2 
     +                        + 0.5*x**4)*(x**2-2.5)
          mu(3,3) = mu(3,3) + weight*voorf*(35./8.-3.5*x**2 
     +                        + 0.5*x**4)**2
        

          if (nv.eq.25600) then
            xpoints(i) = x
            cxi(i) = voorf/dx
          endif

        enddo

        voorf = 1.e19*den(is,ic)*m(is)* 8. / (3.*sqrt(pi))
        mu(2,1) = mu(1,2)
        mu(3,1) = mu(1,3)
        mu(3,2) = mu(2,3)
        meanmuo = 0.
        do k = 1, 3
          do l = 1, 3
            mu(k,l) = mu(k,l) * voorf
            meanmuo = meanmuo + abs(mu(k,l))
          enddo
        enddo
        meanmuo = meanmuo/9.
        nlerr = .false.
        do k = 1, 3
          do l = 1, 3
            if (abs(mu(k,l)-muo(k,l)).gt.abs(epsn*mu(k,l))) then
              if (nv.gt.100) then
                write(*,*) 'WARNING: refinement for viscosity'//
     +                     ' integral needed'
                write(*,'(1x,a8,i1,a1,i1,a4,1p,e15.8)')
     .            '%err.mu(',k,',',l,') = ',
     1           (mu(k,l)-muo(k,l))/mu(k,l)*100.0_R8
                write(*,*) 'voorf =',voorf
                write(*,*) 'nv=',nv
              endif
*    ..neglect errors in minority coefficients
              if (abs(mu(k,l)).gt.(epsn*meanmuo)) then
                nlerr = .true.
              endif
            endif
          enddo
        enddo
        if (nlerr) then
          do k = 1, 3
            do l = 1, 3
              muo(k,l) = mu(k,l)
            enddo
          enddo
          nv = nv*2
          if (nv.gt.4e8) then
            open(16,file='viscos.dat',form='formatted')
            do i = 0, 25600
              write(16,*) xpoints(i), cxi(i)
            enddo
            close(16)
            call perr(4)
          endif
          goto 3000
        endif
         

      return 
      end subroutine viscos


