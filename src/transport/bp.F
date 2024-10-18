      subroutine bp(nsdecl,nzdecl,nsl,ncl,neofrc,t,m,den,zsp,ntau,
     +     xi,la,lab,epsn,nenergy,nleg,nreg,bgradp,fc,fm,mmx,dsi,
     +     epar,coeff)
c--------------------------------------------------------------------
c     input  nsdecl : the maximum number of species
c            nzdecl : the maximum number of charge states
c            nsl     : the number of species
c            ncl     : array(nsdecl+1) the number of charge states
c                      of every species
c            neofrc  : logical if true the matrices are not 
c                      newly calculated 
c            t       : array(nsdecl+1) the temperature of the species
c                      in kev
c            m       : array(nsdecl+1) the mass of the species in kg.
c            den     : array(nsdecl+1,nzdecl) density of every component
c                      in units 10^19 m^-3
c            zsp     : array(nsdecl+1,nzdecl) the charge of every comp.
c                      normalized to e 
c            ntau     : array(nsdecl+1,nsdecl+1) the collision freq. 
c                      m_i n_i / ntau_ij
c            xi      : array(nsdecl+1,nzdecl) the relative weigth of 
c                      every charge state 
c            la      : the friction coefficients
c            lab     : the friction coefficients
c            epsn    : desired accuracy
c            nenergy : integer if 1 energy scattering is used
c                      in the calculation of the viscosity
c            nleg    : the number of Laguerre polynomials 
c            nreg    : a parameter that forces a certain regime
c                      in the calculation of the viscosity
c                      0 weigthed regime valid for all collision
c                        frequencies
c                      1 force the banana regime
c                      other force the Pfirsch-Schlueter regime
c            bgradp  : the flux surface quantity n in grad theta
c                      the inner product of the unit vector in 
c                      the direction of the magnetic field and 
c                      the gradient of the poloidal angle
c            fc      : the number of passing particles
c            fm      : the coefficients necessary to calculate
c                      the viscosity in the Pfirsch-Schlueter 
c                      regime 
c            mmx     : the number of fm coefficients
c            dsi     : the thermodynamic force, for a definition
c                      see routine neoart
c     output coeff   : normalized transport coefficients
c
c     the routine calls the following subroutines
c     perr    : for error handling
c     viscos  : to calculate the viscosity
c     ludcmp  : to do the lu decomposition
c     lubksb  : to solve the matrix equation
c--------------------------------------------------------------------
      use b2mod_types
      use b2mod_constants
      implicit none

      integer :: nsdecl,nzdecl

      integer ::
     *     nsl, ncl(0:nsdecl), i, j, mt, nreg, 
     +         indx(0:3*nsdecl+2), nleg, mmx, k, l
      real (kind=R8) ::
     *       t(0:nsdecl), m(0:nsdecl), den(0:nsdecl,0:nzdecl-1),
     +       ntau(0:nsdecl,0:nsdecl), aaa(3,3,0:nsdecl,0:nsdecl), 
     +       d, akp(3,3), la(3,3,0:nsdecl), 
     +       lab(3,3,0:nsdecl,0:nsdecl), visc(3,3,0:nsdecl,0:nzdecl-1),
     +       xi(0:nsdecl,0:nzdecl-1), fc, fm(mmx), bgradp,
     +       zsp(0:nsdecl,0:nzdecl-1), epar,
     +       ainv(0:3*nsdecl+2,0:3*nsdecl+2), 
     +       ds(0:nsdecl,0:nzdecl-1,2), coeff(0:nsdecl,0:nzdecl-1,4), 
     +       bt(0:3*nsdecl+2), upar(0:nsdecl,0:nzdecl-1,3),
     +       dsi(0:nsdecl,0:nzdecl-1,2)
      logical neofrc


      integer :: nenergy
      real (kind=r8) :: epsn, mu(3,3)

c     define the proper thermodynamic forces for the calculation below
      do i = 0, nsl-1
        do j = 0, ncl(i)-1
          do k = 1, 2
            ds(i,j,k) = t(i)*dsi(i,j,k)/zsp(i,j)
          enddo
        enddo
      enddo

c     check whether the matrices have to be recalculated
      if (.not.neofrc) then 

c     calculate the viscosity
      do i = 0, nsl-1
        do j = 0, ncl(i)-1
          do k = 1, 3
            do l = 1, 3
              visc(k,l,i,j) = 0.0_R8
            enddo
          enddo
          call viscos(nsdecl, nzdecl, nsl, i, j, nreg, bgradp, fc, 
     +                fm, mmx, nenergy, epsn, ntau, m, t, xi, den, 
     +                mu)
          do k = 1, nleg
            do l = 1, nleg 
              visc(k,l,i,j) = mu(k,l)
            enddo
          enddo
        enddo
      enddo

c     calculate the a matrices
      do i = 0, nsl-1 
        do j = 0, ncl(i)-1
          do k = 1, 3
            do l = 1, 3
              aaa(k,l,i,j) = 0.0_R8
            enddo
          enddo
          do k = 1, 3
            aaa(k,k,i,j) = 1.0_R8
          enddo
          do k = 1, nleg
            do l = 1, nleg
              aaa(k,l,i,j) = visc(k,l,i,j) - xi(i,j)*la(k,l,i)
            enddo
          enddo
c         do the lu decomposition
          call ludcmp(aaa(:,:,i,j),2,2,indx,d)
          do k = 1, 3
            do l = 1, 3
              akp(k,l) = 0.0_R8
            enddo
          enddo
          do k = 1, 3
            akp(k,k) = 1.0_R8
          enddo
          do k = 1, 3
            call lubksb(aaa(:,:,i,j),2,2,indx,akp(1,k))
          enddo
          do k = 1, 3
            do l = 1, 3
              aaa(k,l,i,j) = akp(k,l)
            enddo
          enddo

        enddo
      enddo
 
c     now calculate the matrix for inversion
      do i = 0, 3*nsl-1
        do j = 0, 3*nsl-1
          ainv(i,j) = 0.0_R8
        enddo
      enddo
      do i = 0, 3*nsl-1
        ainv(i,i) = 1.0_R8
      enddo
      do i = 0, nsl-1
        do k = 1, 3
          do l = 1, 3
            akp(k,l) = 0.0_R8
          enddo
        enddo
        do k = 1, 3
          do l = 1, 3
            do j = 0, ncl(i)-1
              akp(k,l) = akp(k,l) + xi(i,j)**2*aaa(k,l,i,j)
            enddo
          enddo
        enddo
        do j = 0, nsl-1
          do k = 1, nleg
            do l = 1, nleg
              do mt = 1, 3
                ainv(i*3+k-1,j*3+l-1) = ainv(i*3+k-1,j*3+l-1)
     +             - akp(k,mt)*lab(mt,l,i,j)
              enddo
            enddo
          enddo
        enddo
      enddo

c     calculate the lu decomposition
      call ludcmp(ainv,3*nsl-1,3*nsdecl+2,indx,d)  

c     end the if statement that determines whether the matrices are 
c     recalculated
      endif

c     now calculate the right hand side    
      do i = 0, nsl-1
        do k = 1, 3 
          bt(i*3+k-1) = 0.0_R8
          if (.not.(k.gt.nleg)) then
            do j = 0, ncl(i)-1
              bt(i*3+k-1) = bt(i*3+k-1)+(ev*1.0e19_R8)*xi(i,j)*zsp(i,j)*
     +                        den(i,j)*aaa(k,1,i,j)*epar
            enddo
            do l = 1, 3
              do mt = 1, 2
                do j = 0, ncl(i)-1
                  bt(i*3+k-1) = bt(i*3+k-1) - xi(i,j)*
     +              aaa(k,l,i,j)*visc(l,mt,i,j)*ds(i,j,mt)
                enddo
              enddo
            enddo
          endif
        enddo
      enddo


      call lubksb(ainv,3*nsl-1,3*nsdecl+2,indx,bt)

c     now calculate the indivudual particle velocities
      do i = 0, nsl-1
        do k = 1, 3 
          akp(k,1) = 0.0_R8
          if (.not.(k.gt.nleg)) then
            do l = 1, 3
              do j = 0, nsl-1
                akp(k,1) = akp(k,1)+lab(k,l,i,j)*bt(j*3+l-1)
              enddo
            enddo
          endif
        enddo
        do j = 0, ncl(i)-1
          do k = 1, 3
            upar(i,j,k) = 0.0_R8
            if (.not.(k.gt.nleg)) then
              do l = 1, 3
                akp(l,2) = 0.0_R8
                do mt = 1, 2
                  akp(l,2) = akp(l,2)+visc(l,mt,i,j)*ds(i,j,mt)
                enddo
                upar(i,j,k) = upar(i,j,k) + aaa(k,l,i,j)*(
     +            xi(i,j)*akp(l,1)-akp(l,2))
              enddo
              upar(i,j,k) = upar(i,j,k) + aaa(k,1,i,j)*
     +            (ev*1.0e19_R8)*zsp(i,j)*den(i,j)*epar
            endif
          enddo
        enddo
      enddo


c     finally calculate the coefficients coeff
      do i = 0, nsl-1
        do k = 0, ncl(i)-1

          coeff(i,k,3) = upar(i,k,1)
          do mt = 1, 3
            if (mt.ne.3) 
     +        upar(i,k,mt) = upar(i,k,mt) + ds(i,k,mt)
            coeff(i,k,1) = coeff(i,k,1) + visc(1,mt,i,k)*
     +                     upar(i,k,mt)
            coeff(i,k,2) = coeff(i,k,2) + visc(2,mt,i,k)*
     +                     upar(i,k,mt)
          enddo
          coeff(i,k,4) = upar(i,k,1)

        enddo
      enddo
 
      return
      end subroutine bp
