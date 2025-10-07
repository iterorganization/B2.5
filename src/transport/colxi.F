
      subroutine colxi(nsdecl,nzdecl,nsl,ncl,zsp,den,t,m,ntau,xi)
c--------------------------------------------------------------------
c     this subroutine calculates the collision frequencies and 
c     the weighting factors xi
c
c     input    nsdecl  : maximum number of species
c              nzdecl  : maximum number of charge states
c              nsl     : actual number of species
c              ncl     : array(nsdecl+1) gives the actual amount of 
c                        charge states per species
c              zsp    : array(nsdecl+1,nzdecl) the charge number
c              den    : array(nsdecl+1,nzdecl) the density in 10^19 m^-3
c              t      : array(nsdecl+1) temperature of species
c              m      : array(nsdecl+1) mass of species
c     output   ntau    : array(nsdecl+1,nsdecl+1) collision frequency 
c                       weighted over charge states.
c              xi     : relative weight of every charge state
c--------------------------------------------------------------------
      
      use b2mod_types
      implicit none

      integer ::
     *     nsdecl,nzdecl,i,j,nsl,ncl(0:nsdecl),k,l
      real (kind=R8) ::
     *        den(0:nsdecl,0:nzdecl-1), t(0:nsdecl), m(0:nsdecl),
     &        dene, lnab, ntau(0:nsdecl,0:nsdecl),
     &        xi(0:nsdecl,0:nzdecl-1), tcons, 
     &        zsp(0:nsdecl,0:nzdecl-1)
 
c     calculate the collision times
      tcons = 8.74202e6_R8
c     the collision time below is the double average (over charge 
c     states) of the density divided by the collision time (equation 
c     a3 and a4 of Houlberg, Phys. Plasmas 4 3230 (1997))
c     first determine the electron density
      dene = 0.0_R8
      do i = 0, nsl-1
        if (zsp(i,0).lt.0.0_R8) dene = den(i,0)
      enddo
      if (dene.eq.0.0_R8) then
        do i = 0, nsl-1
          do j = 0, ncl(i)-1
            dene = dene + den(i,j)*zsp(i,j)
          enddo
        enddo
      endif

      do i = 0, nsl-1
        do j = 0, nsl-1

c         calculate the Coulomb logarithm 
          if ((zsp(i,0).lt.0.0_R8).and.(zsp(j,0).lt.0.0_R8)) then
c           electron-electron
            lnab = 14.9_R8 - 0.5_R8*log(0.1_R8*dene)+log(t(i))
          else
            if (zsp(i,0).lt.0.0_R8) then
c             electron-ion
              lnab = 15.2_R8-0.5_R8*log(0.1_R8*dene) + log(t(i))
            else
              if (zsp(j,0).lt.0.0_R8) then
c               ion-electron
                lnab = 15.2_R8-0.5_R8*log(0.1_R8*dene) + log(t(j))
              else
c               ion-ion
                lnab = 17.3_R8-0.5_R8*log(0.1_R8*dene)+
     +                 1.5_R8*log(0.5_R8*(t(i)+t(j)))
              endif
            endif
          endif

          ntau(i,j) = 0.0_R8
          do k = 0, ncl(i)-1
            do l = 0, ncl(j)-1
              ntau(i,j) = ntau(i,j) + den(i,k)*den(j,l)*zsp(i,k)**2
     +                   *zsp(j,l)**2
            enddo
          enddo

          ntau(i,j) = ntau(i,j)*tcons*lnab*sqrt(m(i)/t(i)**3)
 
        enddo
      enddo
   
c     calculate the xi coefficients
      do i = 0, nsl-1
        xi(i,ncl(i)-1) = 0.0_R8
        do j = 0, ncl(i)-1
          xi(i,ncl(i)-1) = xi(i,ncl(i)-1) + den(i,j)*zsp(i,j)**2
        enddo
        do j = 0, ncl(i)-1
          xi(i,j) = den(i,j)*zsp(i,j)**2 / xi(i,ncl(i)-1)
        enddo
      enddo

      return
      end subroutine colxi
