
      subroutine ps(nsdecl,nzdecl,nsl,ncl,neofrc,ntau,xi,la,lab,mas,
     +              temp,den,zsp,sigma,uaioo,rnq,uai,coeff)
c--------------------------------------------------------------------
c     the routine that calculates the Pfirsch-Schlueter contribution.
c
c     input   nsdecl    : maximum number of species
c             nzdecl    : maximum number of charge states
c             nsl     : number of species
c             ncl     : array(nsdecl+1) number of charge states per species
c             neofrc : logical if true then the friction matrices
c                      are not newly calculated
c             ntau    : array(nsdecl+1,nsdecl+1) normalized collision frequen-
c                      cies (m_a n_a / ntau_ab)
c             xi     : array(nsdecl+1,nzdecl) relative weight of the charge
c                      states
c             la     : friction coefficients
c             lab    : friction coefficients
c             mas    : array(nsdecl+1) mass of the species (kg)
c             temp   : array(nsdecl+1) temperature in keV
c             den    : array(nsdecl+1,nzdecl) density in 10^19 m^-3
c             zsp    : charge normalized to e, array(nsdecl+1,nzdecl) 
c             sigma  : control parameters that determine which of 
c                      the coupling terms are taken into account. 
c                      for all coupling terms all sigmas are 1
c                      see routine neoart for definition
c             uaioo  : the thermodynamic forces. array(nsdecl+1,nzdecl,3)
c                      see routine neoart for definition
c             rnq    : length of the field line. in circular geometry
c                      this can be approximated by major radius 
c                      times safety factor 
c     output  uai    : the velocities
c             coeff  : normalized transport coefficients 
c
c     the routine calls the following subroutines
c     perr   : error handling
c     penq   : to calculate the coupling coefficients
c     ludcmp : to calculate the lu decomposition
c     lubksb : to calculate the solution of the matrix equation   
c--------------------------------------------------------------------
      use b2mod_types
      use b2mod_constants
     . ,only: qe
      implicit none

      integer ::
     *     nsdecl, nzdecl, nsl, indx(0:2*nsdecl+1), m, n, ncl(0:nsdecl), 
     +     mt, nt, kt, lt, ms, i,j,k,l
      real (kind=R8) ::
     *       pab(2,2,0:nsdecl,0:nsdecl),qab(2,2,0:nsdecl,0:nsdecl),
     +       la(3,3,0:nsdecl),lab(3,3,0:nsdecl,0:nsdecl),
     +       ntau(0:nsdecl,0:nsdecl),xi(0:nsdecl,0:nzdecl-1),
     +       kab(2,2,0:nsdecl,0:nsdecl),oa(2,2,0:nsdecl),
     +       ha(2,3,0:nsdecl),gk(0:nsdecl),
     +       aaa(0:2*nsdecl+1,0:2*nsdecl+1),a1,a2,a3,lg,da,dnom,
     +       rm(2,2,0:nsdecl,0:nzdecl-1),rmo(2,2,0:nsdecl),
     +       rmoi(2,2,0:nsdecl),dv(2,3,0:nsdecl,0:nsdecl),
     +       ev(2,3,0:nsdecl,0:nsdecl),fv(2,3,0:nsdecl,0:nsdecl),
     +       d,b(0:2*nsdecl+1),uo(0:nsdecl,3),uoo(0:nsdecl,3),
     +       sor(2,0:nsdecl,0:nzdecl-1),soro(2,0:nsdecl),
     +       soroo(2,0:nsdecl),uaio(0:nsdecl,0:nzdecl-1,3),
     +       uai(0:nsdecl,0:nzdecl-1,3),sigma(4),mas(0:nsdecl),
     +       temp(0:nsdecl),den(0:nsdecl,0:nzdecl-1),
     +       uaioo(0:nsdecl,0:nzdecl-1,2),dum1,dum2,pp(2,2),
     +       qq(2,2),ta,tb,ma,mb,lga,
     +       akeep(0:2*nsdecl+1,0:2*nsdecl+1),bkeep(0:2*nsdecl+1),
     +       sol(0:2*nsdecl+1),btus(0:2*nsdecl+1),f1,f2,rnq,
     +       coeff(0:nsdecl,0:nzdecl-1,4),usolo(0:nsdecl,3),
     +       rmoo(2,2,0:nsdecl),rmooo(2,2,0:nsdecl),
     +       zsp(0:nsdecl,0:nzdecl-1)
      logical nliter, nlchan, neofrc      
      nliter = .false.
      nlchan = .false.
            
c     some coefficients that appear in the equations
      a1 = 4.0_R8/25.0_R8
      a2 = 16.0_R8/175.0_R8
      a3 = (8.0_R8/35.0_R8)**2

c     since the thermodynamic forces are proportional to t/ z 
c     the uaio quantities are multiplied by this number.
      do i = 0, nsl-1
        do j = 0, ncl(i)-1 
          do k = 1, 2
            uaio(i,j,k) = temp(i)*uaioo(i,j,k) /(zsp(i,j))
          enddo
          uaio(i,j,3) = 0.0_R8
        enddo
      enddo

c     if the routine is called with the same densities and temperatures
c     then do not repeat the calculation of the matrix
      if (.not.neofrc) then

c     the coefficient that determines the influence of the field
c     line length
      lg = -1./(rnq)**2

c     calculated the friction coefficients
      do i = 0, nsl-1
        do  j = 0, nsl-1
          ta = temp(i)
          tb = temp(j)
          ma = mas(i)
          mb = mas(j)       
          call penq(ta,tb,ma,mb,pp,qq)
          do k = 1, 2
            do l = 1, 2
              pab(k,l,i,j) = pp(k,l)
              qab(k,l,i,j) = qq(k,l)
            enddo
          enddo
        enddo
      enddo

c     calculate the kab and oa matrices
      do i = 0, nsl-1
        dum1 = 1./(temp(i))
        oa(1,1,i) = 0.0_R8
        oa(1,2,i) = 0.0_R8
        oa(2,1,i) = 0.0_R8
        oa(2,2,i) = 0.0_R8
        do j = 0, nsl-1
          dum2 = 1./(temp(j))
          kab(1,1,i,j) = -a1*dum1*sigma(3)*ntau(i,j)*
     +                   pab(1,1,i,j)
          kab(1,2,i,j) = - a2*dum2*sigma(1)*ntau(i,j)*
     +                   qab(1,2,i,j)
          kab(2,1,i,j) = a2*dum1*sigma(2)*ntau(i,j)*
     +                   (pab(2,1,i,j)+pab(1,1,i,j))
          kab(2,2,i,j) = a3*dum2*ntau(i,j)*sigma(4)*
     +                   (qab(2,2,i,j)+sigma(1)*qab(1,2,i,j))
          oa(1,1,i) = oa(1,1,i) + a1*dum1*sigma(3)*
     +                ntau(i,j)*pab(1,1,i,j)
          oa(1,2,i) = oa(1,2,i) - a2*dum1*sigma(1)*
     +                ntau(i,j)*pab(1,2,i,j)
          oa(2,1,i) = oa(2,1,i) - a2*dum1*sigma(2)*
     +                ntau(i,j)*(pab(2,1,i,j) +
     +                pab(1,1,i,j))
          oa(2,2,i) = oa(2,2,i) + a3*dum1*ntau(i,j)*sigma(4)*
     +                (pab(2,2,i,j)+ sigma(1)*pab(1,2,i,j))
        enddo
      enddo

c     now calculate the ha coefficients
      do i = 0, nsl-1
        do k = 1, 2
          do l = 1, 3
            ha(k,l,i) = 0.0_R8
            do j = 1, 2
              ha(k,l,i) = ha(k,l,i) - oa(k,j,i)*
     +                    la(j+1,l,i)
            enddo
          enddo
        enddo
      enddo

      do i = 0, nsl-1
        da = ha(1,2,i)*ha(2,3,i)-ha(1,3,i)*ha(2,2,i)
        ta = ha(1,2,i)+ha(2,3,i)
        lga = lg*mas(i)*(qe*1000.0_R8*1.0e19_R8**2)
        do j = 0, ncl(i)-1
          dnom = (xi(i,j)/den(i,j))/ (lga**2 + (xi(i,j)/
     +       (den(i,j)))**2*lga*ta+(xi(i,j)/(den(i,j)))**4*da)
          dum1 = (xi(i,j)/den(i,j))**2
          rm(1,1,i,j) =  dnom*(lga+dum1*ha(2,3,i))
          rm(1,2,i,j) = -dnom*dum1*ha(1,3,i)
          rm(2,1,i,j) = -dnom*dum1*ha(2,2,i)
          rm(2,2,i,j) =  dnom*(lga+dum1*ha(1,2,i))
        enddo
      enddo
 
      do i = 0, nsl-1
        gk(i) = 0.0_R8
        do k = 1, 2
          do l = 1, 2
            rmo(k,l,i) = 0.0_R8
            rmoo(k,l,i) = 0.0_R8
            rmooo(k,l,i) = 0.0_R8
          enddo
        enddo
        do j = 0, ncl(i)-1
          dum1 = xi(i,j)**2/den(i,j)
          dum2 = xi(i,j)*dum1/den(i,j)
          gk(i) = gk(i) + dum1
          do k = 1, 2
            do l = 1, 2
              rmo(k,l,i) = rmo(k,l,i)+ xi(i,j)*rm(k,l,i,j)
              rmoo(k,l,i) = rmoo(k,l,i)+dum1*rm(k,l,i,j)
              rmooo(k,l,i) = rmooo(k,l,i)+dum2*rm(k,l,i,j)
            enddo
          enddo
        enddo
      enddo
 
      do i = 0, nsl-1
        dnom = 1.0_R8/(rmo(1,1,i)*rmo(2,2,i)-
     +                 rmo(1,2,i)*rmo(2,1,i))
        rmoi(1,1,i) =  dnom*rmo(2,2,i)
        rmoi(1,2,i) = -dnom*rmo(1,2,i)
        rmoi(2,1,i) = -dnom*rmo(2,1,i)
        rmoi(2,2,i) =  dnom*rmo(1,1,i)
      enddo

      do i = 0, nsl-1
        do j = 0, nsl-1

          do k = 1, 2
            do l = 1, 3
              dv(k,l,i,j) = 0.0_R8
              ev(k,l,i,j) = 0.0_R8
              fv(k,l,i,j) = 0.0_R8
              do m = 1, 2
                dv(k,l,i,j) = dv(k,l,i,j) + 
     +                        kab(k,m,i,j)*la(m+1,l,j)
                ev(k,l,i,j) = ev(k,l,i,j) + 
     +                        oa(k,m,i)*lab(m+1,l,i,j)
                do n = 0, nsl-1
                  fv(k,l,i,j) = fv(k,l,i,j) + 
     +               gk(n)*kab(k,m,i,n)*lab(m+1,l,n,j)
                enddo
              enddo
            enddo
          enddo

        enddo
      enddo

c     build the matrix
      do i = 0, nsl-1
        do j = 0, nsl-1
          do k = 1, 2
            do l = 1, 2
              dum1 = 0.0_R8
              dum2 = 0.0_R8
              do n = 1, 2
                dum1 = dum1 + rmo(k,n,i)*fv(n,l+1,i,j)
     +                 +rmoo(k,n,i)*ev(n,l+1,i,j)
                do mt = 1, 2
                  do nt = 1, 2
                     dum2 = dum2 + rmo(k,n,i)*dv(n,mt+1,i,j)*
     +                 rmoo(mt,nt,j)*rmoi(nt,l,j)
                  enddo
                enddo
              enddo
              aaa(2*i+k-1,2*j+l-1) = - dum1 - dum2

              dum1 = 0.0_R8
              do n = 0, nsl-1
                do kt = 1, 2
                  do lt = 1, 2 
                    do mt = 1, 2
                      dum2 = 0.0_R8
                      do nt = 1, 2
                        do ms = 1, 2
                          dum2 = dum2 + rmoo(lt,nt,n)*
     +                      rmoi(nt,ms,n)*rmoo(ms,mt,n)
                        enddo
                      enddo
                      dum2 = rmooo(lt,mt,n)-dum2
                      dum1 = dum1 + rmo(k,kt,i)*dv(kt,lt+1,i,n)*
     +                    dum2*ev(mt,l+1,n,j)
                    enddo
                  enddo
                enddo
              enddo
              aaa(2*i+k-1,2*j+l-1) = aaa(2*i+k-1,2*j+l-1)
     +           - dum1
            enddo
          enddo
        enddo
      enddo

      do i = 0, 2*nsl-1
        aaa(i,i) = aaa(i,i) + 1.0_R8
      enddo

      if (nlchan) then
      endif

      if (nliter) then
        do i = 0, 2*nsl-1
          do j = 0, 2*nsl-1
            akeep(i,j) = aaa(i,j)
          enddo
        enddo
      endif

c     now do the lu decomposition 
      call ludcmp(aaa,2*nsl-1,2*nsdecl+1,indx,d)

c     end the if statement that determines whether the same 
c     matrix is to be used. 
      endif

      do i = 0, nsl-1
        do k = 1, 3
          uo(i,k) = 0.0_R8
          uoo(i,k) = 0.0_R8
          do n = 0, ncl(i)-1
            uo(i,k) = uo(i,k) + xi(i,n)*uaio(i,n,k)
            uoo(i,k)=uoo(i,k)+xi(i,n)**2*uaio(i,n,k)/
     +         (den(i,n))
          enddo
        enddo
      enddo

c     calculate the right hand side solutions
      do i = 0, nsl-1
        do j = 0, ncl(i)-1
          do k = 1, 2
            sor(k,i,j) = 0.0_R8
            do l = 1, 3
              do mt = 1, 2
                sor(k,i,j) = sor(k,i,j) - (xi(i,j)/den(i,j))*
     +            rm(k,mt,i,j)*ha(mt,l,i)*uaio(i,j,l)
                do n = 0, nsl-1
                  sor(k,i,j) = sor(k,i,j) + rm(k,mt,i,j)*
     +            (fv(mt,l,i,n)*uo(n,l)+dv(mt,l,i,n)*uoo(n,l)+
     +            (xi(i,j)/den(i,j))*ev(mt,l,i,n)*uo(n,l))
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      do i= 0, nsl-1
        do k = 1, 2
          soro(k,i) = 0.0_R8
          soroo(k,i) = 0.0_R8
          do j = 0, ncl(i)-1
            soro(k,i) = soro(k,i) + xi(i,j)*sor(k,i,j)
            soroo(k,i) = soroo(k,i) + xi(i,j)**2*sor(k,i,j)/
     +                   den(i,j)
          enddo
        enddo
      enddo

      do i = 0, nsl-1
        do k = 1, 2
          b(2*i+k-1) = soro(k,i)
          do n = 0, nsl-1
            do mt = 1, 2
              do nt = 1, 2
                b(2*i+k-1) = b(2*i+k-1) + rmo(k,mt,i)*
     +            dv(mt,nt+1,i,n)*soroo(nt,n)
                do kt = 1, 2
                  do lt = 1, 2
                    b(2*i+k-1) = b(2*i+k-1) - rmo(k,mt,i)*
     +                dv(mt,nt+1,i,n)*rmoo(nt,kt,n)*
     +                rmoi(kt,lt,n)*soro(lt,n)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
           

      if (nlchan) then
      endif

      if (nliter) then
        do i = 0, 2*nsl-1
          bkeep(i) = b(i)
        enddo
      endif
      
c     calculate the solution through backsubstition   
      call lubksb(aaa,2*nsl-1,2*nsdecl+1,indx,b)
 
      if (nliter) then
        do i = 0, 2*nsl-1
          sol(i) = b(i)
        enddo

        do i = 0, 2*nsl-1
          btus(i) = bkeep(i)
          do j = 0, 2*nsl-1
            btus(i) = btus(i)-akeep(i,j)*sol(j)
          enddo
        enddo
        call lubksb(aaa,2*nsl-1,2*nsdecl+1,indx,btus)
        do i = 0, 2*nsl-1
          sol(i) = sol(i)+btus(i)
        enddo
        do i = 0, 2*nsl-1
          b(i) = sol(i)
        enddo
      endif

c     now calculate the velocities
      do i = 0, nsl-1
        do j = 0, ncl(i)-1
          do k = 1, 2
            uai(i,j,k+1) = sor(k,i,j)
            do mt = 1, 2
              do nt = 1, 2
                uai(i,j,k+1) = uai(i,j,k+1) + rm(k,mt,i,j)*
     +             rmoi(mt,nt,i) * (b(2*i+nt-1) - soro(nt,i))
              enddo
            enddo
            do kt = 1, 2
              do lt = 1, 2
                dum1 = 0.0_R8
                if (kt.eq.lt) dum1 = dum1 + xi(i,j)/den(i,j)
                do mt = 1, 2
                  dum1 = dum1 - rmoi(kt,mt,i)*rmoo(mt,lt,i)
                enddo
                do mt = 1, 2 
                  do n = 0, nsl-1
                    uai(i,j,k+1) = uai(i,j,k+1) + rm(k,kt,i,j)*
     +                dum1*ev(lt,mt+1,i,n)*b(2*n+mt-1)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
                  
  
      do i = 0, nsl-1
        do j = 0, ncl(i)-1
          uai(i,j,1) = uaio(i,j,1)
          uai(i,j,2) = uai(i,j,2)+uaio(i,j,2)
          usolo(i,1) = uo(i,1)
          usolo(i,2) = uo(i,2) + b(2*i)
          usolo(i,3) = b(2*i+1)
        enddo
      enddo
 
c     calculate the friction 
      do i = 0, nsl-1
        f1 = 0.0_R8
        f2 = 0.0_R8
        do j = 0, nsl-1
          do k = 1, 3
            f1 = f1 + lab(1,k,i,j)*usolo(j,k)
            f2 = f2 + lab(2,k,i,j)*usolo(j,k)
          enddo
        enddo
        do j = 0, ncl(i)-1
          coeff(i,j,1) = 0.0_R8
          coeff(i,j,2) = 0.0_R8
          do k = 1, 3
            coeff(i,j,1) = coeff(i,j,1) + xi(i,j)*(la(1,k,i)*
     +                     uai(i,j,k))
            coeff(i,j,2) = coeff(i,j,2) + xi(i,j)*(la(2,k,i)*
     +                     uai(i,j,k))
          enddo
          coeff(i,j,1) = coeff(i,j,1) + xi(i,j)*f1
          coeff(i,j,2) = coeff(i,j,2) + xi(i,j)*f2
        enddo
      enddo
 
      return
      end subroutine ps
